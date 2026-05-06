#!/bin/bash

usage(){
	echo "This script creates a coverage prfile, in the form of a BedGraph file, for each BAM file in the directory and subdirectory of the script."
	echo "It automatically creates a track containig the type, name and description for each BAM file"
	echo ""
	echo "USAGE"
	echo "./coverage_script.sh [OPTIONS]"
	echo ""
	echo "OPTIONS"
	echo "-h          Prints the help for this script"
	exit 0
}

while getopts ":h" opt; do # checks which possible options are beeing selected
	case $opt in
		h) # check if the help option is present
			usage # calls the function with the instructions
			;;
		\?) # checks if an invalid option was passed
			echo "Error: Invalid option -$OPTARG" >&2; usage ;;
		:) # checks if all options requiring an argument actually recived one
			echo "Error: Option -$OPTARG requires an argument" >&2; usage ;;
	esac
done

shift $((OPTIND - 1)) # makes sure the options are not counted for positional arguments

parent_folders=$(find | egrep "\.bam$" | sed -r "s/[^/]*$//" | sort | uniq) # gets the parent paths the BAM files

for parent_folder in $parent_folders; do
	BAM_files=$(ls $parent_folder | egrep \.bam$) # names of the BAM files in the current folder
	BAM_paths=$(find $parent_folder | egrep \.bam$) # paths of the BAM files in the current folder
	parent_folder_name=$(echo $parent_folder | egrep -o [^/]+/$ | sed "s/\///") # gets the name of the current folder

	for BAM in $BAM_files; do # cycles through all BAM file names in the current directory
		current_BAM_path="${parent_folder}${BAM}" # single BAM path of one of the BAM in the directory
		track_desc="'coverage track for ${BAM}'" # variable containing the description to show in the UCSC genome browser
		output_bg_path=$(echo "${parent_folder}${BAM}_cov.bg" | sed "s/.bam//") # path of the bed graph file

		# creation of the coverage profile for each BAM file
		bedtools genomecov -ibam $current_BAM_path -bg -trackopts "name=${BAM} description=${track_desc}" -max 100 > $output_bg_path
	done
done
