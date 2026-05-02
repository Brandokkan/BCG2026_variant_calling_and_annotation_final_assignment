#!/bin/bash

usage(){
	echo "This script creates quality control of all the BAM file in the directory and sub directories of the script."
	echo "Qualimap is used for the coverage, fastqc for the quality of the reads and multiqc is used to create a final report"
	echo "This script automatically handles the necessary arguments and compatibility needed for the tools used for quality control"
	echo "By default, the whole genome is checked for QC"
	echo ""
	echo "USAGE"
	echo "./multiQC_script.sh [options] <qualimap_tsv_name>"
	echo ""
	echo "qualimap_tsv_name:   The name of the tsv file used by qualimap to recognize the BAM files. There must be one per directory containing BAM files."
	echo "                     The tsv file name follows this pattern: (name of the paret directory)_(input name)"
	echo "                     To be recognized, the tsv file must follow the pattern. If no tsv file is present, one is created with this pattern and"
	echo "                     points to all the BAM files in the directory."
	echo ""
	echo "OPTIONS"
	echo "-h                        Prints the help for this script"
	echo "-b <target_region.BED>    Tells the script to check only certain target regions for QC. the target regions must be specified in a BED file"
	echo "-i <multiqc_report_name>  Adds a specific name to the final multiqc report"
	exit 0
}

BED_present=false # used to check if the BED file containing the regions of interest for qualimap is present
multiQC_name=multiQC_report # name of the final multiQC report
# $1 name of tsv file containig the BAM files information
while getopts ":hb:i:" opt; do # checks which possible options are beeing selected
	case $opt in
		h) # check if the help option is present
			usage # calls the function with the instructions
			;;
		b)
			BED_present=true
			BED_file=$OPTARG # name (or pattern) of the BED file containing the regions of interest
			;;
		i)
			multiQC_name=$OPTARG
			;;
		\?) # checks if an invalid option was passed
			echo "Error: Invalid option -$OPTARG" >&2; usage ;;
		:) # checks if all options requiring an argument actually recived one
			echo "Error: Option -$OPTARG requires an argument" >&2; usage ;;
	esac
done

shift $((OPTIND - 1)) # makes sure the options are not counted for positional arguments


parent_folders=$(find | egrep \.bam$ | sed -r "s/[^/]*$//" | sort | uniq) # gets the parent paths the BAM files

for parent_folder in $parent_folders; do
	BAM_files=$(ls $parent_folder | egrep \.bam$) # names of the BAM files in the current folder
	BAM_paths=$(find $parent_folder | egrep \.bam$) # paths of the BAM files in the current folder
	parent_folder_name=$(echo $parent_folder | egrep -o [^/]+/$ | sed "s/\///") # gets the name of the current folder

	# checks if the tsv file exists. If not, it creates it with the BAM files in the current folder
	if [[ ! -e "${parent_folder}${parent_folder_name}_${1}" ]]; then
		for BAM_file in $BAM_files; do
			echo -e "${BAM_file}\t${parent_folder}${BAM_file}\t${parent_folder_name}" >> "${parent_folder}${parent_folder_name}_${1}"
		done
	fi

	# cheks if the user selected the option to use a file containig the regions of interest and then chooses the correct option
	if [[ $BED_present == true ]]; then
		qualimap multi-bamqc -d "${parent_folder}${parent_folder_name}_${1}" -gff $BED_file -outdir "${parent_folder}qualimap_files" -r
	else
		qualimap multi-bamqc -d "${parent_folder}${parent_folder_name}_${1}" -outdir "${parent_folder}qualimap_files" -r
	fi

	# for multiqcc compatibility, the names of the stats folder must be changed
	stats_folders=$(find $parent_folder | egrep _stats$) # paths of the folders created by qualimap for the single BAM file reports
	for stats_folder in $stats_folders; do
		new_name=$(echo $stats_folder | sed "s/_stats//") # new name without the final _stat part
		mv $stats_folder $new_name # renaming
	done

	mkdir "${parent_folder}fastq_files"
	fastqc -o "${parent_folder}fastq_files" $BAM_paths

done

multiqc -i $multiQC_name -x nohup.out . # creates a multiQC report combining togheter the outputs of the previous QC analisys
