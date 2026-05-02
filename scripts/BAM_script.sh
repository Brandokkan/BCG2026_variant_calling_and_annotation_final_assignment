#!/bin/bash

# script that creates BAM and BAM index files from compressed paired-end FastQ files

usage(){ # help function
	echo "./BAM_script.sh [options] ref_name pattern"
	echo ""
	echo "This script creates BAM and BAM index files from compressed pair-end FastQ files."
	echo "It looks in the directory and sub-directories of the script for pair-end fq.gz fiels"
	echo "To work, the pair-end files must end with 1.fq.gz and 2.fq.gz respectly"
	echo ""
	echo "Positional arguments"
	echo "ref_name: the common name of the bowtie2 reference genome used"
	echo "pattern: tells the script to take into consideration only the fq.gz files that contain this pattern in their directory"
	echo ""
	echo "Options"
	echo "-h: prints the help function for this script"
	echo "-i reference_fasta: tells the script that the bowtie2 reference genome is not already present and must be built first."
	echo "                    needs the name of the reference fasta file"
	exit 0
}

index_present=true # used to check if a bowtie index of the reference genome is already present
# $1 è il nome dell'index bowtie, $2 è il nome comune di TUTTI file fq.gz da analizzare
while getopts ":hi:" opt; do # checks which possible options are beeing selected
	case $opt in
		h) # check if the help option is present
			usage # calls the function with the instructions
			;;
		i) # checks if the index option is present
			index_present=false
			reference_fa=$OPTARG # sets the name of the reference genome fasta file
			;;
		\?) # checks if an invalid option was passed
			echo "Error: Invalid option -$OPTARG" >&2; usage ;;
		:) # checks if all options requiring an argument actually recived one
			echo "Error: Option -$OPTARG requires an argument" >&2; usage ;;
	esac
done

shift $((OPTIND - 1)) # makes sure the options are not counted for positional arguments


if [ $index_present == false ]; then
	bowtie2-build $reference_fa $1 # builds the bowtie2 index for the reference genome choosen
fi

# variable with the paths of all pair-end compressed fq files (one name for each pair)
fq_paths=$(find | egrep \.fq\.gz$ | sed -r "s/(1|2)\.fq\.gz$//" | sort | uniq -d)

for path_name in $fq_paths; do
	cur_path=$(echo $path_name | sed -r "s/[^/]*$//") # extracts the path without the name of the fq file
	parent_folder=$(echo $cur_path | egrep -o [^/]+/$ | sed "s/\///") # gets the name of the folder containing the fq file

	# checks which family member is the BAM
	if [[ $path_name == *"421"* ]]; then
		BAM_name="${parent_folder}_child"
		id="S421"
	elif [[ $path_name == *"422"* ]]; then
		BAM_name="${parent_folder}_father"
		id="S422"
	else
		BAM_name="${parent_folder}_mother"
		id="S423"
	fi

	if [[ $path_name == *$2* ]]; then # checks if the fq file contain the specified pattern adn then creates the BAM and BAM index
		final_BAM_path="${cur_path}${BAM_name}.bam"
		bowtie2 -p 10 -x $1 -1 "${path_name}1.fq.gz" -2 "${path_name}2.fq.gz" --rg-id $id --rg "SM:${BAM_name}" | samtools view -b | samtools sort -o $final_BAM_path
		samtools index $final_BAM_path
	fi
done
