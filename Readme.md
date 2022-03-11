# eggd_vep

## What does this app do?

Annotates a vcf using [Variant Effect Predictor](https://github.com/Ensembl/ensembl-vep). Default docker image used v104.3.

## What are typical use cases for this app?
This app was designed to annotate vcfs with specified fields based on provided annotation.


## What data are required for this app to run?
- An input vcf to annotate
- VEP docker image (`vep_docker`)
- VEP plugins (`vep_plugin_config`) - VEP version specific:
    - plugin_config.txt
- VEP reference files (`vep_refs`):
    - Homo_sapiens.GRCh37.dna.toplevel.fa.gz
    - Homo_sapiens.GRCh37.dna.toplevel.fa.gz.fai
    - Homo_sapiens.GRCh37.dna.toplevel.fa.gz.gzi
    - homo_sapiens_refseq_vep_105_GRCh37.tar.gz
- Annotation configuration file (`config_file`):
    - json file providing information about annotations and plugins.
    - Example config file:
  	```
    {
        "custom_annotations": [
        {
            "name": "ClinVar",
            "type": "vcf",
            "annotation_type": "exact",
            "force_coordinates": "0",
            "vcf_fields": "CLNSIG,CLNREVSTAT,CLNDN",
            "resource_files": [
                {
                "file_id": "project-Fkb6Gkj433GVVvj73J7x8KbV:file-G6BY8X0433GbVBG06pFPvjp7",
                "index_id": "project-Fkb6Gkj433GVVvj73J7x8KbV:file-G6BY8Pj433GjPF7j6kK5QZ75"
                }
            ]
        }
    ],
        "plugins": [
            {
            "name": "CADD",
            "pm_file": "project-Fkb6Gkj433GVVvj73J7x8KbV:file-G6BYQGQ433Gg0KfZF5xPfv1X",
            "resource_files": [
                {
                "file_id": "project-Fkb6Gkj433GVVvj73J7x8KbV:file-G5KjGZ84YV3xFfv221vJQjQZ",
                "index_id": "project-Fkb6Gkj433GVVvj73J7x8KbV:file-G5KjKf84YV3y39z97vbXJ60b"
                },
                {
                "file_id": "project-Fkb6Gkj433GVVvj73J7x8KbV:file-G6BYB5Q433GvQGgF58976p7J",
                "index_id": "project-Fkb6Gkj433GVVvj73J7x8KbV:file-G6BYB8j433GyJ3f13z8GB4ZV"
                }
            ]
        }
    ]
    }

	```
> All the annotation sources above are specific for GRCh37 and are set us default for the app - this app can run with the equivalent annotation sources for GRCh38.

__This app uses the follow tools which are app assets:__
* htslib (v1.14)
* bedtools (v2.30.0)

## What are the optional inputs for this app?
- buffer_size [default = 500] : to allow for parallelisation the app recognised the instance type and splits annotation in the amount of available cores. The buffer size is the amount of variants VEP will annotate per core.

> For larger vcfs please consider about the appropriate  instance type and buffer size to use.
## What does this app output?
- Annotated vcf with the specified fields mentioned above.

## Notes
- Designed to be used as part of Helios workflow for processing Solid Cancer data for GRCh37 so all defaults are based on that.
- This app uses a buffer_size of 500 variants and parallelised the maximum number of cores available.As a default this app runs using mem1_ssd2_v2_x4 which translates to 4 cores. This was chosen with the solid cancer vcfs in mind as mentioned in the inputs section the buffer size can be given as an argument and the instance changed at runtime.
