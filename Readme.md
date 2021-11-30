# eggd_vep

## What does this app do?

Annotates a vcf using [Variant Effect Predictor](https://github.com/Ensembl/ensembl-vep). Default docker image used v104.3.

## What are typical use cases for this app?
This app was designed to annotate vcfs with specified fields. 

Annotate against specified refseq transcripts with
- gene symbol
- variant class
- variant consequence
- exon number
- HGVS c.
- HGVS p.
- gnomAD AF
- gnomADg AF (genomes)
- CADD PHRED
- dbSNP
- ClinVar
- ClinVar - Clinical Indication
- ClinVar - Clinical Significance
- Transcript Feature
   

## What data are required for this app to run?
- An input vcf to annotate
- VEP docker image (`vep_docker`)
- VEP plugins (`vep_plugins`)
    - CADD.pm
    - plugin_config.txt
- VEP reference files (`vep_refs`):
    - Homo_sapiens.GRCh37.dna.toplevel.fa.gz
    - Homo_sapiens.GRCh37.dna.toplevel.fa.gz.fai
    - Homo_sapiens.GRCh37.dna.toplevel.fa.gz.gzi
    - homo_sapiens_refseq_vep_104_GRCh37.tar.gz
* VEP annotation sources (`vep_annotation`):
   - Cosmic Coding Variants VCF (v94)
      - CosmicCodingMuts_v94_grch37.normal.vcf.gz
      - CosmicCodingMuts_v94_grch37.normal.vcf.gz.tbi
  - Cosmic NonCoding Variants VCF (v94)
      - CosmicNonCodingVariants_v94_grch37.normal.vcf.gz
      - CosmicNonCodingVariants_v94_grch37.normal.vcf.gz.tbi
  - ClinVar VCF (2021113)
      - clinvar_20211113_withChr.vcf.gz
      - clinvar_20211113_withChr.vcf.gz.tbi
  - CADD (v1.6)
      - gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz
      - gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz.tbi
      - cadd_whole_genome_SNVs_GRCh37.tsv.gz
      - cadd_whole_genome_SNVs_GRCh37.tsv.gz.tbi
  
  
> All the annotation sources above are specific for GRCh37 and are set us default for the app - this app can run with the equivalent annotation sources for GRCh38.

__This app uses the follow tools which are app assets:__
* htslib (v1.14)
* bedtools (v2.30.0)
## What does this app output?
- Annotated vcf with the specified fields mentioned above.

## Limitations
- Designed to be used as part of Helios workflow for processing Solid Cancer data for GRCh37 so all defaults are based on that.
- This app uses a buffer_size of 500 variants and parallelised using 4 cores. This was chosen with the solid cancer vcfs in mind and may need to be tweaked or changed for bigger vcfs.
  
