# BCG2026 genomics final assignment
This repository contains the scripts and data used for the final assignment
of Professor Chiara for the 2026 BCG genomics and transcriptomics course.

## Data
Contains the final output data used for the analysis.

### QC report
Contains the multiQC report produced by multiQC_script.sh, which is used for quality control of the BAM files.

### UCSC genome browser images
Images representing the coverage tracks of trios with a disease and the corresponding variant in the UCSC genome browser.

### Unclear case
VEP and web-based VEP of trio 4 produced by a variation of the pathogenic_va_script.sh to include a variant less rare than the threshold.
See the final project report for further details.

### VEP command line files
Final VCF produced by pathogenic_va_script.sh using the command line of VEP present on the UNIMI server.

### VEP web files
txt versions of the VCF files in "VEP command line files". Further annotated using the web version of VEP.
NOTE: The trio_1 and trio_4 VEP web versions are not present because no variants passed all pathogenic_va_script.sh filters.

### bedGraph files
Bedgraph files containing the read coverage profile for all the family members in the trios.

## Scripts
Contain the scripts used for the analysis. The scripts can take -h as an argument to print the help function.

### BAM_script.sh
Creates BAM files from paired-end FASTQ files.

### coverage_script.sh
Creates a bedgraph file containing the coverage profile for each BAM file.

### multiQC_script.sh
Creates a multiQC report of all the BAM files, merging the QC from qualimap and fastqc.

### pathogenic_va_script.sh
Creates a final VCF file containing all the VEP-annotated variants that passed all the filters to be considered pathogenic.
It also outputs intermediate files, most of which are not present in this repo.

## Author
Brando S.V.A
