# eggd_vep

## What does this app do?

Annotates a vcf using [Variant Effect Predictor](https://github.com/Ensembl/ensembl-vep). Default docker image used v105.0.

## What are typical use cases for this app?
This app was designed to annotate vcfs with specified fields based on provided annotation.

A variable level of annotation can be achieved by different combinations of custom annotation and vep plugins, in addition to the required VEP cache annotation bundle.

## What data are required for this app to run?
- An input vcf to be annotated (`vcf`)
- Annotation configuration file (`config_file`):
    - json file providing information about annotations and plugins.
    - Example config file:
  	```
    {
            "config_information":{
                "genome_build": "GrCh37",
                "assay":"TWE",
                "config_version": "1.0.0"
            },

            "vep_resources":{
                "vep_docker":"project-Fkb6Gkj433GVVvj73J7x8KbV:file-G8V2Vz0433Gp5bYPF2f6vg9X",
                "vep_cache":"project-Fkb6Gkj433GVVvj73J7x8KbV:file-G8V4bGj433Gz96K3Fb1VfbG3",
                "plugin_config":"project-Fkb6Gkj433GVVvj73J7x8KbV:file-G8V57Yj433Gfg3vF9jPq1ZFk",
                "reference_fasta":"project-Fkb6Gkj433GVVvj73J7x8KbV:file-G6BYyyj4YV3pYBkgFVGP2K4P",
                "reference_fai":"project-Fkb6Gkj433GVVvj73J7x8KbV:file-G6P32kj4YV3y58KyP4k4qG2p",
                "reference_gzi":"project-Fkb6Gkj433GVVvj73J7x8KbV:file-G6BYz104YV3X5qp463K5b5vp",
                "ref_bcftools":"project-Fkb6Gkj433GVVvj73J7x8KbV:file-F3zxG0Q4fXX9YFjP1v5jK9jf"

            },
            "custom_annotations": [
                {
                    "name": "ClinVar",
                    "type": "vcf",
                    "annotation_type": "exact",
                    "force_coordinates": "0",
                    "vcf_fields": "CLNSIG,CLNREVSTAT,CLNDN",
                    "required_fields":"ClinVar,ClinVar_CLNSIG,ClinVar_CLNDN",
                    "resource_files": [
                        {
                        "file_id":"project-Fkb6Gkj433GVVvj73J7x8KbV:file-G6BY8X0433GbVBG06pFPvjp7",
                        "index_id":"project-Fkb6Gkj433GVVvj73J7x8KbV:file-G6BY8Pj433GjPF7j6kK5QZ75"
                        }
                    ]
                }
            ],
            "plugins": [
                {
                    "name": "SpliceAI",
                    "pm_file": "project-G86K7XQ4jKXPgK4Z8Zj585Zj:file-G8YXQQQ4jKX4bVz53zK7VJX5",
                    "required_fields":"SpliceAI_pred_DS_AG,SpliceAI_pred_DS_AL,SpliceAI_pred_DS_DG,SpliceAI_pred_DS_DL",
                    "suffix": "cutoff=0.5",
                    "resource_files": [
                    {
                        "file_id": "project-Fkb6Gkj433GVVvj73J7x8KbV:file-G9FFfBj433GV57Zf8ZvxfbBg",
                        "index_id": "project-Fkb6Gkj433GVVvj73J7x8KbV:file-G9FFfz0433GZ1X13FqjjQJFF",
                        "prefix": "snv="
                    },
                    {
                        "file_id": "project-Fkb6Gkj433GVVvj73J7x8KbV:file-G9FF6zj433Gv5jxkFqpK6p7J",
                        "index_id": "project-Fkb6Gkj433GVVvj73J7x8KbV:file-G9FF6qQ433GxvVp497P5V05X",
                        "prefix": "indel="
                    }
                ]
            }
        ]
    }


	```
> In theory the app should run with any combination of annotation, please bear in mind the instance used if passing big datasets


## What are the optional inputs for this app?
- Amount of variants VEP will annotate per core(`buffer_size`) [default = 500] : to allow for parallelisation the app recognises the instance type and splits annotation in the amount of available cores.
- A panel bed file to filter the vcf on (`panel_bed`)
- A list of transcripts to filter on (`transcript_list`). One transcript per line. VEP annotates with all possible transcripts and if this list is passed it filters on the given transcript list.
- A boolean flag of whether to normalise the input vcf or not (`normalise`) [ default=true ].

__This app uses the following tools which are app assets:__
* htslib (v1.14)
* bedtools (v2.30.0)



> For larger vcfs please consider the appropriate instance type and buffer size to use.
## What does this app output?
- Annotated (and if requested, filtered) vcf.

## Notes
- This app uses a buffer_size of 500 variants and parallelised the maximum number of cores available. As a default, this app runs using mem1_ssd1_v2_x16 which translates to 16 cores. This was chosen to speed up set up.
- The default behaviour of this app is to normalise the input vcf as all default annotation used is also normalised to be able to appropriately compare and annotate the vcf. The `normalise` option was built to ensure compatibility with copy-number vcfs which do not require normalisation.
- If there are additional flags, the `additional_flags` section in the config file should contain the flag name and the `additional_fields` section in the config file should contain the output field name for that flag. Example flags and output fields are linked [here](https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html)
