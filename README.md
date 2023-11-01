# NestedProkkaFix
Script to check for falsely split orthologues in roary outputs

# Usage
````
Nested_prokka_annotation_fix.py [-h] -o OUTPUT [-p PRESENCE_PERCEN] -r REFERENCE (-i ROARY_CSV | -id ROARY_WITH_DETAILS) (-gff INPUT_GFF_FOLDER | -dict GFF_DICTIONARY)
````
## Parameters Description



options:
- `h, --help`  :          show this help message and exit
- `o OUTPUT, --output OUTPUT`:
                        Output direcotry. (default: None)
- `p PRESENCE_PERCEN, --presence_percen PRESENCE_PERCEN`:
                        Percentage of isolates a gene must be present. (default: 0.95)
- `r REFERENCE, --reference REFERENCE`:
                        Reference accession. (default: None)
- `i ROARY_CSV, --roary_csv ROARY_CSV`:
                        Roary presence and absence matrix. (default: None)
- `id ROARY_WITH_DETAILS, --roary_with_details ROARY_WITH_DETAILS`:
                        Roary presence and absence matrix with annotationd details added. (default: None)
- `gff INPUT_GFF_FOLDER, --input_gff_folder INPUT_GFF_FOLDER`:
                        path to folder containing all the prokka gff file (default: None)
- `dict GFF_DICTIONARY, --gff_dictionary GFF_DICTIONARY`:
                        Preprocessed .JSON of GFF annotations. (default: None)