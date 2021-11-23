#!/bin/bash
#
# Performs annotation and filtering of given mutect2 VCF to produce multiple annotated VCFs
# Performs filtering and annotation of given pindel VCF to produce annotated VCF
# n.b. filtering for pindel is a combination of a bed file of regions, removal of all 1-2bp
# insetions and filtering by transcript with vep
# Generates an excel workbook to aid variant interpretation
# Also generates a workaround for BSVI mis-handling multiallelics

function annotate_vep_vcf {
	# Function to run VEP for annotation on given VCF file

	# Inputs:
	# $1 -> input vcf 
	# $2 -> name for output vcf

	input_vcf="$1"
	output_vcf="$2"
	
	# fields to filter on
	# hard coded in function for now, can be made an input but all are the same
	filter_fields="SYMBOL,VARIANT_CLASS,Consequence,EXON,HGVSc,HGVSp,gnomAD_AF,gnomADg_AF,CADD_PHRED,Existing_variation,ClinVar,ClinVar_CLNDN,ClinVar_CLNSIG,COSMIC,Feature"

	# find clinvar vcf, remove leading ./
	clinvar_vcf=$(find ./ -name "clinvar_*.vcf.gz" | sed s'/.\///')

	# find cosmic coding and non-coding vcfs
	cosmic_coding=$(find ./ -name "CosmicCodingMuts*.vcf.gz" | sed s'/.\///')
	cosmic_non_coding=$(find ./ -name "CosmicNonCodingVariants*.vcf.gz" | sed s'/.\///')

	# find CADD files, remove leading ./
	cadd_snv=$(find ./ -name "*SNVs.tsv.gz")
	cadd_indel=$(find ./ -name "*indel.tsv.gz")

	# find gnomad files, remove leading ./
  	gnomad_genome_vcf=$(find ./ -name "gnomad.genomes.*.vcf.gz" | sed s'/.\///')

	time docker run -v /home/dnanexus:/opt/vep/.vep \
	ensemblorg/ensembl-vep:release_104.3 \
	./vep -i /opt/vep/.vep/"${input_vcf}" -o /opt/vep/.vep/"${output_vcf}" \
	--vcf --cache --refseq --exclude_predicted --symbol --hgvs --af_gnomad \
	--check_existing --variant_class --numbers \
	--offline \
	--custom /opt/vep/.vep/"${clinvar_vcf}",ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
	--custom /opt/vep/.vep/"${cosmic_coding}",COSMIC,vcf,exact,0,ID \
	--custom /opt/vep/.vep/"${cosmic_non_coding}",COSMIC,vcf,exact,0,ID \
  --custom /opt/vep/.vep/"${gnomad_genome_vcf}",gnomADg,vcf,exact,0,AF,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH \
	--plugin CADD,/opt/vep/.vep/"${cadd_snv}",/opt/vep/.vep/"${cadd_indel}" \
	--fields "$filter_fields" \
	--no_stats
}

main() {
	set -e -x -v -o pipefail

	mark-section "downloading inputs"
	time dx-download-all-inputs --parallel

	# array inputs end up in subdirectories (i.e. ~/in/array-input/0/), flatten to parent dir
	find ~/in/vep_plugins -type f -name "*" -print0 | xargs -0 -I {} mv {} ~/in/vep_plugins
	find ~/in/vep_refs -type f -name "*" -print0 | xargs -0 -I {} mv {} ~/in/vep_refs
	find ~/in/vep_annotation -type f -name "*" -print0 | xargs -0 -I {} mv {} ~/in/vep_annotation

	# move annotation sources to home
	mv ~/in/vep_annotation/* /home/dnanexus/

	mark-section "annotating"
	
	# vep needs permissions to write to /home/dnanexus
	chmod a+rwx /home/dnanexus
	
	# extract vep reference annotation tarball to /home/dnanexus
	time tar xf /home/dnanexus/in/vep_refs/*.tar.gz -C /home/dnanexus

	# place fasta and indexes for VEP in the annotation folder
	mv /home/dnanexus/in/vep_refs/*fa.gz* ~/homo_sapiens_refseq/104_GRCh37/

	# place plugins into plugins folder
	mkdir ~/Plugins
	mv ~/in/vep_plugins/* ~/Plugins/


	# load vep docker
	docker load -i "$vep_docker_path"

	#annotate
	output_vcf="${vcf_prefix}_annotated.vcf"
	print ($output_vcf)
	annotate_vep_vcf "$vcf" "$output_vcf"

	mark-section "uploading output"

  # upload output files

	# Upload output vcf
    annotated_filtered_vcf=$(dx upload $output_vcf --brief)

    # The following line(s) use the utility dx-jobutil-add-output to format and
    # add output variables to your job's output as appropriate for the output
    # class.  Run "dx-jobutil-add-output -h" for more information on what it
    # does.

    dx-jobutil-add-output annotated_filtered_vcf "$annotated_filtered_vcf" --class=file
	
	
	mark-success
}