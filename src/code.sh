#!/bin/bash
#
# Annotates input vcf with variant effect predictor (VEP) outputing the filter fields as seen in line 17

_annotate_vep_vcf () {
	# Function to run VEP for annotation on given VCF file

	# Inputs:
	# $1 -> input vcf
	# $2 -> name for output vcf

	input_vcf="$1"
	output_vcf="$2"

	# fields to annotate with.

	fields="SYMBOL,VARIANT_CLASS,Consequence,EXON,HGVSc,HGVSp,gnomAD_AF,Existing_variation,Feature"
	for entry in $(jq -r '.custom_annotations,.plugins| .[].required_fields' "$config_file_path" );
	do
		fields+=','$entry;
	done

	# Get number of cores/threads to use in --fork option
	echo "Number of forks: $FORKS"
	echo "Buffer size used: $buffer_size"

	# --exclude_null_allelels is used with --check-existing to prevent multiple COSMIC id's
	# being added to the same variant.
	# the buffer size is chosen based on the average size of the input VCF

	/usr/bin/time -v docker run -v /home/dnanexus:/opt/vep/.vep \
	${VEP_IMAGE_ID} \
	./vep -i /opt/vep/.vep/"${input_vcf}" -o /opt/vep/.vep/"${output_vcf}" \
	--vcf --cache --refseq --exclude_predicted --symbol --hgvs --af_gnomad \
	--check_existing --variant_class --numbers \
	--offline --exclude_null_alleles \
	$ANNOTATION_STRING $PLUGIN_STRING --fields "$fields" \
	 --buffer_size "$buffer_size" --fork "$FORKS" \
	--no_stats --compress_output bgzip --shift_3prime 1
}

_filter_vep_vcf () {
	# Function to filter annotated VCF with VEP to retain variants with AF < 0.10 in gnomAD, for
	# a gene symbol being present and against a given list of transcripts

	# Inputs:
	# 	$1 -> input vcf (should be output vcf of annotation)
	# 	$2 -> name for output_vcf
	#	$3 -> file containing one transcript per line

	input_vcf="$1"
	output_vcf="$2"
	transcripts="$3"

	# Create string filter with NM ids an no version
	# VEP expects Feature match <transcript_id> to match partially
	# Separate transcript separated by an "or"
	transcript_list=$(for tr in $(less $transcripts);do echo -n "Feature match $tr\. or ";done)

	/usr/bin/time -v docker run -v /home/dnanexus:/opt/vep/.vep \
	${VEP_IMAGE_ID}  \
	./filter_vep -i /opt/vep/.vep/"$input_vcf" \
	-o /opt/vep/.vep/"$output_vcf" --only_matched --filter \
	"${transcript_list%????}"
}

_format_annotation () {
    # Formats the annotation part of the command given a config file

	# Inputs:
	# $1 -> input config file
    local file=$1

    ANNOTATION_STRING=""

    for annotation in $(jq -c '.custom_annotations[]' "$file")
    do
        ANNOTATION_STRING+=" --custom "

        # Get the file name using the file id and add to the command list
        dx_file_id=$(jq -r '.resource_files[0].file_id' <<< "$annotation")
        dx_name=$(dx describe "$dx_file_id" --json | jq -r '.name')
        ANNOTATION_STRING+="/opt/vep/.vep/${dx_name},"

        # Adds required vep arguments
        ANNOTATION_STRING+=$(jq -j  '[
            .name,.type,.annotation_type,.force_coordinates,(.vcf_fields // empty)
        ] | map(tostring) | join(",")' <<< "$annotation")
    done

}


_format_plugins () {
    # Formats the plugin part of the command given a config file

	# Inputs:
	# $1 -> input config file
    local file=$1

    PLUGIN_STRING=""

    for plugin in $(jq -c '.plugins[]' "$file")
    do
        # Add the required flag and name for VEP plugins
        local plugin_name
        plugin_name=$(jq -r '.name' <<< "$plugin")

        PLUGIN_STRING+=" --plugin ${plugin_name},"

        for plugin_file in $(jq -rc '.resource_files[]' <<< "$plugin");
        do
            local prefix
            prefix=$(jq -j '(.prefix // empty)' <<< "$plugin_file")

            # Get the file name using the file id
            # Add it the command list
            dx_file_id=$(jq -r '.file_id' <<< "$plugin_file")
            dx_name=$(dx describe "$dx_file_id" --json | jq -r '.name')
            PLUGIN_STRING+="${prefix}/opt/vep/.vep/${dx_name},"
        done

        # Check if there's additional options to append
        if [ "$(jq -r 'has("suffix")' <<< $plugin)" = true ]
        then
            PLUGIN_STRING+=$(jq -r '.suffix' <<< $plugin)
        else
            # Remove trailing comma
            PLUGIN_STRING="${PLUGIN_STRING%?}"
        fi

    done
}


main() {
	set -e -x -v -o pipefail

	FORKS=$(grep -c ^processor /proc/cpuinfo)

	mark-section "downloading inputs"
	time dx-download-all-inputs --parallel

	# Print Config information to logs:
	jq -r '.config_information | to_entries | .[] | .key + ": " + (.value )' "$config_file_path"
	# Download vep resources
	vep_refs=$(jq -r '.vep_resources | .[]' "$config_file_path")
	xargs -P"$FORKS" -n1 dx download <<< $vep_refs

	# Makes the names of the files accessible by their key name in the config
	# Creates variables storing the name of vep_resources
	# i.e. vep_docker=vep_v105.0.tar.gz which makes the name accessible further down
	eval $(jq -r '.vep_resources
		| to_entries | .[] 
		| .key + "=$(dx describe --json " + (.value )+"
		| jq -r '.name')" ' "$config_file_path")


	# Download annotation files
	echo $(date +%T)

	vep_files=$(jq -r '.custom_annotations,.plugins | .[].resource_files[] | .file_id,  (.index_id // empty) ' "$config_file_path")
	xargs -P"$FORKS" -n1 dx download <<< $vep_files

	echo $(date +%T)

	# Download plugin .pm files
	if [[ ! -z $(jq '.plugins[].pm_file' "$config_file_path") ]];
	then
		dx download $(jq -r '.plugins[].pm_file' "$config_file_path")
			# place plugins into plugins folder
		mkdir ~/Plugins
		mv ~/$plugin_config ~/Plugins/
		mv ~/*.pm ~/Plugins/
	else
		echo "No plugin files found."
	fi

	mark-section "pre-annotation filtering & normalisation"

	# Unpack fasta reference
	#fasta_tar=$(dx describe $(jq -r '.vep_resources | .' vep_testing_config_v0.0.json) --json | jq -r '.name')
	tar xzf $ref_bcftools

	# Filter by panel if provided and normalise
	if [ "$panel_bed" ];
	then
		bedtools intersect -header -a "$vcf_path" -b "$panel_bed_path" \
			| bcftools norm -f genome.fa -m -any --keep-sum AD - \
			-o "${vcf_prefix}_normalised.vcf"
	else
		bcftools norm -f genome.fa -m -any --keep-sum AD "$vcf_path"  \
		-o "${vcf_prefix}_normalised.vcf"
	fi

	# If hard filters are passed in the config apply them.
	filter_command=$(jq -r '.quality_filters' "$config_file_path")

	# Test if quality_filters variable is missing or empty in the config
	if test -z "$filter_command" || [ $filter_command == "null" ]
	then
		echo "No filtering commands passed."
		mv "${vcf_prefix}_normalised.vcf" "${vcf_prefix}_temp.vcf"
	else
		echo "Filter commands passed: $filter_command"
		eval $(echo ${filter_command} "${vcf_prefix}_normalised.vcf" -o "${vcf_prefix}_temp.vcf")
	fi

	mark-section "annotating"

	# vep needs permissions to write to /home/dnanexus
	chmod a+rwx /home/dnanexus

	# extract vep reference annotation tarball to /home/dnanexus
	time tar xf $vep_cache -C /home/dnanexus

	# Place fasta and indexes for VEP in the annotation folder
	chmod a+rwx /home/dnanexus/*fa.gz*

	# Find the annotation cache folder dynamically to allow different versions of VEP and genome to be used
	# This will only work with refseq annotation caches currently as this is passed in the VEP command
	# We can consider changing this if non-refseq files are required but it needs the --refseq flag to become optional
	cache_path=$(find homo_sapiens_refseq -mindepth 1 -maxdepth 1 -type d -name "*_GRCh3*" )
	echo "$cache_path"
	mv /home/dnanexus/*fa.gz* ~/"${cache_path}"

	# load vep docker
	docker load -i "$vep_docker"
	VEP_IMAGE_ID=$(docker images --format="{{.Repository}} {{.ID}}" | grep "^ensemblorg" | cut -d' ' -f2)

	# Create annotation and plugin strings
	_format_annotation "$config_file_path"
	_format_plugins "$config_file_path"

	# Annotate
	annotated_vcf="${vcf_prefix}_temp_annotated.vcf.gz"
	echo $annotated_vcf
	_annotate_vep_vcf "${vcf_prefix}_temp.vcf" "$annotated_vcf"

	# Filter vcf by chosen transcript
	output_vcf="${vcf_prefix}_annotated.vcf"

	if [ "$panel_bed" ];
	then
		# Extract the transcripts noted in the panel bed and remove their version
		cut -f 4 $panel_bed_path | sort | uniq |cut -d '.' -f 1 > transcripts.tsv
		# Filter annotated output with chosen transcript
		_filter_vep_vcf "${vcf_prefix}_temp_annotated.vcf.gz" ${output_vcf} transcripts.tsv
		bgzip ${output_vcf}
	elif [ "$transcript_list" ];
	then
		# Ensure transcript list has no versions
		sort "$transcript_list_path" | uniq | cut -d '.' -f 1 > transcripts.tsv
		# Filter annotated output with chosen transcript
		_filter_vep_vcf "${vcf_prefix}_temp_annotated.vcf.gz" ${output_vcf} transcripts.tsv
		bgzip ${output_vcf}
	else
		mv "${vcf_prefix}_temp_annotated.vcf.gz" "${output_vcf}.gz"
	fi

	# Upload output vcf
	mark-section "uploading output"
	annotated_vcf=$(dx upload "${output_vcf}.gz" --brief)
	dx-jobutil-add-output annotated_vcf "$annotated_vcf" --class=file

	mark-success
}