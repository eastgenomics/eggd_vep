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
	# hard coded in function for now, can be made an input but all are the same
	#filter_fields="SYMBOL,VARIANT_CLASS,Consequence,EXON,HGVSc,HGVSp,gnomAD_AF,gnomADg_AF,CADD_PHRED,Existing_variation,ClinVar,ClinVar_CLNDN,ClinVar_CLNSIG,COSMIC,Feature"

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
	$ANNOTATION_STRING $PLUGIN_STRING \
	 --buffer_size "$buffer_size" --fork "$FORKS" \
	--no_stats --compress_output --shift_3prime
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

	# array inputs end up in subdirectories (i.e. ~/in/array-input/0/), flatten to parent dir
	find ~/in/vep_plugins -type f -name "*" -print0 | xargs -0 -I {} mv {} ~/in/vep_plugins
	find ~/in/vep_refs -type f -name "*" -print0 | xargs -0 -I {} mv {} ~/in/vep_refs

	# move annotation sources and input vcf to home
	mv ~/in/vcf/* /home/dnanexus/

	# Download annotations
	echo $(date +%T)

	vep_files=$(jq -r '.[] | .[].resource_files[] | .file_id,  (.index_id // empty) ' "$config_file_path")
	xargs -P"$FORKS" -n1 dx download <<< $vep_files

	echo $(date +%T)

	# Download plugin .pm files
	dx download $(jq -r '.plugins[].pm_file' "$config_file_path")

	mark-section "annotating"

	# vep needs permissions to write to /home/dnanexus
	chmod a+rwx /home/dnanexus

	# extract vep reference annotation tarball to /home/dnanexus
	time tar xf /home/dnanexus/in/vep_refs/*.tar.gz -C /home/dnanexus

	# place fasta and indexes for VEP in the annotation folder
	chmod a+rwx /home/dnanexus/in/vep_refs/*fa.gz*

	# Find the annotation cache folder dynamically to allow different versions of VEP and genome to be used
	# This will only work with refseq annotation caches currently as this is passed in the VEP command
	# We can consider changing this if non-refseq files are required but it needs the --refseq flag to become optional
	cache_path=$(find homo_sapiens_refseq -mindepth 1 -maxdepth 1 -type d -name "*_GRCh3*" )
	echo "$cache_path"
	mv /home/dnanexus/in/vep_refs/*fa.gz* ~/"${cache_path}"

	# place plugins into plugins folder
	mkdir ~/Plugins
	mv ~/in/vep_plugin_config/* ~/Plugins/
	mv ~/*.pm  ~/Plugins/

	# load vep docker
	docker load -i "$vep_docker_path"
	VEP_IMAGE_ID=$(docker images --format="{{.Repository}} {{.ID}}" | grep "^ensemblorg" | cut -d' ' -f2)

	# Create annotation and plugin strings
	_format_annotation "$config_file_path"
	_format_plugins "$config_file_path"

	# Annotate
	output_vcf="${vcf_prefix}_annotated.vcf.gz"
	echo $output_vcf
	_annotate_vep_vcf "$vcf_name" "$output_vcf"


	# Upload output vcf
	mark-section "uploading output"
	annotated_vcf=$(dx upload $output_vcf --brief)
	dx-jobutil-add-output annotated_vcf "$annotated_vcf" --class=file

	mark-success
}