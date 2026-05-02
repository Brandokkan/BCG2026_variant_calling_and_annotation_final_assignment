#!/bin/bash

usage(){
	echo "This script creates a vcf file containing all candidates pathogenic variants starting from the BAM files of a family sample."
	echo "It also creates intermediate files like the complete vcf, the inheritance compatible vcf and the VEP format of the former for each family sample."
	echo "Variant calling is done using freebayes, while variant annotation is done using VEP."
	echo "The above files are produced for each sub-directory of the directory in which the script is located containing at least one BAM file."
	echo ""
	echo "USAGE"
	echo "./pathogenic_va_script.sh [OPTIONS] <reference_genome> <mode_of_inheritance_file>"
	echo ""
	echo "reference_genome:             Reference genome fasta file used for variant calling."
	echo ""
	echo "mode_of_inheritance_file:     tsv file containing the type of hineritance (AR, AD_denovo, AD_inherited) for each family sample."
	echo "                              To work, the tsv must contain the exact name of the folder containing the family BAM files and have on the same line"
	echo "                              the type of inheritence of that family written as in the example above."
	echo "                              for AD_inherited, the affected parent must be specified on the same line (ex. mother_affected). If no parent is"
	echo "                              specified, the script will assume both parents are affected."
	echo ""
	echo "OPTIONS"
	echo "-h                            Prints the help for this script"
	echo ""
	echo "-b <target_region.BED>        Tells the script to consider only variants in specified regions contained in a BED file"
	exit 0
}

# $1 reference genome
# $2 modes of inheritance file
while getopts ":hb:i:" opt; do # checks which possible options are beeing selected
	case $opt in
		h) # check if the help option is present
			usage # calls the function with the instructions
			;;
		b)
			BED_present=true
			BED_file=$OPTARG # name (or pattern) of the BED file containing the regions of interest
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

	vcf_name="${parent_folder_name}.vcf" # name of the file that will contain the variants called with freebayes
	vcf_path="${parent_folder}${vcf_name}"
	freebayes -f $1 -m 20 -C 5 -Q 10 -q 10 --min-coverage 10 $BAM_paths > $vcf_path # variant calling

	bgzip $vcf_path
	vcf_name_comp="${vcf_name}.gz"
	vcf_path_comp="${parent_folder}${vcf_name_comp}"
	bcftools index $vcf_path_comp

	vcf_candidates_name="cand_${vcf_name}" # name of the file that will contain the variants compatible with the inheritance filters
	vcf_candidates_path="${parent_folder}${vcf_candidates_name}"
	vcf_order_file_path="${parent_folder}${parent_folder_name}_vcf_order"
	echo $BAM_files | sort | sed "s/.bam//g" | sed "s/\s/\n/g" > $vcf_order_file_path  # file containing the names of the samples in alphabetical order used for consistent vcf filtering
	type_of_inh_cur_line=$(grep $parent_folder_name $2) # the line in the MOI file that corresponds to the type of inheritance of the curent samples
	if [[ $type_of_inh_cur_line == *AD_denovo* ]]; then # checks what is the type of inheritance and assigns the proper filter variable
		toi_filter='GT[0]!="RR"'
	elif [[ $type_of_inh_cur_line == *AD_inherited* ]]; then
		if [[ $type_of_inh_cur_line == *mother_affected* ]]; then
			toi_filter='GT[0]="RA" && GT[2]!="RR"'
		elif [[ $type_of_inh_cur_line == *father_affected* ]]; then
			toi_filter='GT[0]="RA" && GT[1]!="RR"'
		else
			toi_filter='GT[0]!="RR" && GT[1]!="RR" && GT[2]!="RR"'
		fi
	else
		toi_filter='GT[0]="AA" && GT[1]!="RR" && GT[2]!="RR"'
	fi

	# checks if the user selected the option to only study specific regions and runs the program accordingly
	if [[ $BED_present == true ]]; then
		# filters for variants that match with the selected inheritance filters
		bcftools view -R $BED_file $vcf_path_comp | bcftools view -S $vcf_order_file_path | bcftools view -i "${toi_filter}" | bcftools filter -i 'QUAL>20' -Ov -o $vcf_candidates_path
	else
		bcftools view $vcf_path_comp | bcftools view -S vcf_order | bcftools view -i $toi_filter | bcftools filter -i 'QUAL>20' -Ov -o $vcf_candidates_pat
	fi

	vep_name="${parent_folder_name}.vep" # name of the file that will contain the annotated variants using VEP
	vep_path="${parent_folder}${vep_name}"
	# variant annotation
	vep -i $vcf_candidates_path -o $vep_path --vcf --cache --offline --assembly GRCh37 --dir_cache /data/vep_cache --use_given_ref --mane --pick_allele --af --af_1kg --af_gnomade --max_af --sift b --polyphen b

	vep_filtered_name="${parent_folder_name}_filtered.vep" # name of the file that will contain the candidates pathogenic variants
	vep_filtered_path="${parent_folder}${vep_filtered_name}"
	# filtering for field values that match with pathogenic variants
	filter_vep -i $vep_path -o $vep_filtered_path --filter "(IMPACT is HIGH or (IMPACT is MODERATE and SIFT <= 0.05 and PolyPhen >= 0.15)) and (not MAX_AF or MAX_AF < 0.0001)"
done
