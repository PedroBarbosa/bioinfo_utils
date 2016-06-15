#!/bin/bash

display_usage() { 
printf "	-1st argument must be a file listing all the pairs to process. In this file the Paired end libraries need to be first. If one concatenated reads file is provided in the '.fq' or '.fastq' format the command will be generated too. In this case no matter if you put true/false in the second and third aguments, only assembly command will be generated.
	-2nd argument must be a flag true/false to use paired end reads to generate the command.
	-3rd argument must be a flag true/false to use mate pair reads to generate the command.
	-4th argument must be the number of threads to use.
	-5th argument must be the maximum amount of memory allowed.
	-6th argument must be the size of the k-mer size to start with the assembly. [Even if only the scaffolding is desired, set this argument as you want for this script to work].
	-7th argument is optional. If set to yes, only performs the assemble Available options: [yes|no] Default: [Empty argument, which implies assembly and scaffolding together].
	-8th argument is optional. If set to yes, only performs scaffolding. The prefix needs to be same as the contigs file, normally 'output'. Available options: [yes|no]. Default: [Empty argument, assembly and scaffolding together,if 6th argument not supplied. If 6th supplied as yes, an error will be thrown].
	-9th argument is optional. It refers to the orientation of the mate pair libraries. Available options: [rf|fr]. Default:[rf]
	-10th argument is optional. If set to yes, platanus will use the given mp/pe insert sizes for scaffolding. Available options: [no|yes]. Default: [no]\n\n"
}


function assembly_concatenated_file() {

	command_contigs="$command_contigs -f $1"	

	#check if output file already exists. If so, delete it
        file=command_concatenated_file.sh
        if [ -f "$file" ]; then
                rm $file
        fi

        ##Pass command to script and give permissions to run
        echo -e "#!/bin/bash \n$command_contigs" > ./$file | chmod +x ./$file

}

function assembly_contigs() {
	numb_samples=0
	while read line
	do	
		filename=$line

		#check if reached the MATE PAIR samples, and if not supposed to assemble them, break the loop
		if [[ "$filename" == *"MP"* && "$3" = "true" ]];then
			matepairFlag=true

		elif [[ "$filename" == *"MP"* && "$3" = "false" ]];then
			break
		fi
		
		#add files to the command [even though they are paired end there is no difference for the contig step
		if [ -f "$filename" -a "$matepairFlag" = "false" ]; then
			file=$filename
			let "numb_samples += 1"
			string_samples="$string_samples $file"

		elif [ -f "$filename" -a "$matepairFlag" = "true" ]; then
			file=$filename
			let "numb_samples += 1"
			string_samples="$string_samples $file"
		fi
		
	
	done <"$1"
	numb_samples=$(($numb_samples / 2))
	#add samples to the command
	if [[ "$2" = "true" ]]; then

		command_contigs="$command_contigs $string_samples"
	elif [[ "$2" = "false" ]]; then	
		printf "Are you trying to assemble without paired end reads ? Please set true for the use of this read files\n"
		display_usage
		exit 1
	else
		printf "Please provide valid values for the type of libraries to be used\n"
		display_usage
		exit 1
	fi

}

function scaffolding(){
contigs_file="-c output_contig.fa"
bubbles_file="-b output_contigBubble.fa"
#Paired end and/or MP reads [mp reads are treated as IP because they are already with the correct orientation]. Case letter on IP because it is two different files for each pair
IP=""
OP=""
INS_STRING=""
STDV_STRING=""
numb_samples=0
matepairFlag=false
first_pair=true
second_pair=false
declare -a INSERT_SIZES_PE=('170' '170' '170' '170' '170' '170' '170' '170' '500' '500' '500' '500' '500' '500' '500' '500' '500' '800' '800' '800' '800' '800' '800' '800' '800' '800' '800' '800' '800')
declare -a INSERT_STDV_PE=('30' '30' '30' '30' '30' '30' '30' '30' '60' '60' '60' '60' '60' '60' '60' '60' '60' '90' '90' '90' '90' '90' '90' '90' '90' '90' '90' '90' '90')
declare -a INSERT_SIZES_MP=('2100' '2100' '2100' '2100' '2100' '2100' '2100' '2100' '2100' '2100' '2100' '2100' '5000' '5000' '5000' '5000' '5000' '5000')
declare -a INSERT_STDV_MP=('250' '250' '250' '250' '250' '250' '250' '250' '250' '250' '250' '250' '600' '600' '600' '600' '600' '600')

#mate pair orientatin provided?
if [ -z "$4" ] ; then 
	orientation="outward"
elif [ "$4" = "fr" ]; then
	orientation="inward"
elif [ "$4" = "rf" ]; then
	orientation="outward"
else
	printf "\nPlease provide a correct value for the orientation of the mate pair libraries.\n\n"
	display_usage
	exit 1
fi

while read line
do	
	filename=$line

	#check if reached the MATE PAIR samples, and if not supposed to assemble them, break the loop
	if [[ "$filename" == *"MP"* && "$3" = "true" ]];then
		matepairFlag=true

	elif [[ "$filename" == *"MP"* && "$3" = "false" ]];then
		break
	fi
	
	#paired end
	if [ -f "$filename" -a "$first_pair" = "true" -a "$2" = "true" -a "$matepairFlag" = "false" ]; then
		pair1=$filename
		first_pair=false
		second_pair=true

	elif [ -f "$filename" -a  "$second_pair" = "true" -a "$2" = "true" -a "$matepairFlag" = "false"  ]; then
		pair2=$filename
		first_pair="true"
		second_pair="false"       
		let "numb_samples += 1"	  
		IP="$IP -IP${numb_samples} $pair1 $pair2 "

	elif [ ! -f "$filename" -a "$matepairFlag" = "false"  ]; then
		echo "$filename" is not a file
		continue
	fi
	
	#mate pair
	if [ -f "$filename" -a "$first_pair" = "true" -a "$3" = "true" -a "$matepairFlag" = "true" ]; then
                pair1=$filename
                first_pair=false
                second_pair=true

    elif [ -f "$filename" -a  "$second_pair" = "true" -a "$3" = "true" -a "$matepairFlag" = "true" -a "$orientation" = "inward"  ]; then
                pair2=$filename
                first_pair="true"
                second_pair="false"
                let "numb_samples += 1"
                IP="$IP -IP${numb_samples} $pair1 $pair2 "
	
	elif [ -f "$filename" -a  "$second_pair" = "true" -a "$3" = "true" -a "$matepairFlag" = "true" -a "$orientation" = "outward"  ]; then
                pair2=$filename
                first_pair="true"
                second_pair="false"
                let "numb_samples += 1"
                OP="$OP -OP${numb_samples} $pair1 $pair2 "

    elif [ ! -f "$filename" -a "$matepairFlag" = "false"  ]; then
                echo "$filename" is not a file
                continue

    fi

done < "$1"

#add samples to the command
if [ "$3" = "true" -a "$orientation" = "inward" ]; then
	if [ -z "$5" ];then
		command_scaffolds="$command_scaffolds $contigs_file $bubbles_file $IP"	
	elif [ "$2" = "true" ]; then

	    #add ins and stdv for pe libs
        for (( i = 0; i < ${#INSERT_SIZES_PE[@]}; i++ )); do
            INS_STRING="$INS_STRING -a{$(expr $i + 1) ${INSERT_SIZES_PE[$i]} "
            STDV_STRING="$STDV_STRING -d{$(expr $i + 1) ${INSERT_STDV_PE[$i]} "
        done

        #add ins and stdv for mp libs
        for (( i = 0; i < ${#INSERT_SIZES_MP[@]}; i++ )); do
            INS_STRING="$INS_STRING -a{$(expr ${#INSERT_SIZES_PE[@]} + $i + 1) ${INSERT_SIZES_MP[$i]} "
            STDV_STRING="$STDV_STRING -d{$(expr ${#INSERT_STDV_PE[@]} + $i + 1) ${INSERT_STDV_MP[$i]} "
        done
		command_scaffolds="$command_scaffolds $contigs_file $bubbles_file $IP $INS_STRING $STDV_STRING"

	else
	    #add ins and stdv for mp libs
        for (( i = 0; i < ${#INSERT_SIZES_MP[@]}; i++ )); do
            INS_STRING="$INS_STRING -a{$(expr $i + 1) ${INSERT_SIZES_MP[$i]}"
            STDV_STRING="$STDV_STRING -d{$(expr $i + 1) ${INSERT_STDV_MP[$i]}"
        done
		command_scaffolds="$command_scaffolds $contigs_file $bubbles_file $IP $INS_STRING $STDV_STRING"
	fi

elif [ "$3" = "true" -a "$orientation" = "outward" ]; then
	
	if [ -z "$5" ] ;then
		if [ "$2" = "true" ] ; then
			command_scaffolds="$command_scaffolds $contigs_file $bubbles_file $IP $OP"
		else
			command_scaffolds="$command_scaffolds $contigs_file $bubbles_file $OP"
		fi
	else
		if [ "$2" = "true" ]; then

		    #add ins and stdv for pe libs
            for (( i = 0; i < ${#INSERT_SIZES_PE[@]}; i++ )); do
                INS_STRING="$INS_STRING -a$(expr $i + 1) ${INSERT_SIZES_PE[$i]}"
                STDV_STRING="$STDV_STRING -d$(expr $i + 1) ${INSERT_STDV_PE[$i]}"
            done

            #add ins and stdv for mp libs

            for (( i = 0; i < ${#INSERT_SIZES_MP[@]}; i++ )); do
                INS_STRING="$INS_STRING -a$(expr ${#INSERT_SIZES_PE[@]} + $i + 1) ${INSERT_SIZES_MP[$i]}"
                STDV_STRING="$STDV_STRING -d$(expr ${#INSERT_STDV_PE[@]} + $i + 1) ${INSERT_STDV_MP[$i]}"
            done
			command_scaffolds="$command_scaffolds $contigs_file $bubbles_file $IP $OP $INS_STRING $STDV_STRING"

		else

		    #add ins and stdv for mp libs
            for (( i = 0; i < ${#INSERT_SIZES_MP[@]}; i++ )); do
                INS_STRING="$INS_STRING -a$(expr $i + 1) ${INSERT_SIZES_MP[$i]}"
                STDV_STRING="$STDV_STRING -d$(expr $i + 1) ${INSERT_STDV_MP[$i]}"
            done
			command_scaffolds="$command_scaffolds $contigs_file $bubbles_file $OP $INS_STRING $STDV_STRING"
		fi
		
	fi
	
elif [ "$3" = "false" ]; then	
	printf "\nAre you trying to scaffolding without mate pair reads ? Please set true for the use of this reads files or set the seventh argument to yes, in order to generate the assembly command only. \n\n"
	display_usage
	exit 1

else
	printf "\nPlease provide valid values for the type of libraries to be used.\n\n"
	display_usage
	exit 1
fi

}

function generate_final_file_contigs(){
	#check if output file already exists. If so, delete it
	file=command_${numb_samples}_samples.sh
	if [ -f "$file" ]; then
		rm $file
	fi

	##Pass command to script and give permissions to run
	echo -e "#!/bin/bash \n$command_contigs" > ./$file | chmod +x ./$file
}

function generate_final_file_scaffolds(){
	#check if output file already exists and if it is supposed to append. If so, delete it
	file=command_${numb_samples}_samples.sh
	if [ -f "$file" ] && [ -z $1 ] ; then #$1 refers to the first parameter passed to the function [in this case, if only scaffolding]. This case in supposed to append, because no information was passed in the script arguments
		echo -e "\n$command_scaffolds" >> ./$file
	elif [[ -f "$file" &&  $1 = "true" ]] ; then
		rm $file
		echo -e "#!/bin/bash \n$command_scaffolds" > ./$file | chmod +x ./$file
	else
		echo -e "#!/bin/bash \n$command_scaffolds" > ./$file | chmod +x ./$file
	fi	
}




#########################################################
#################General variables#######################
#########################################################
#exec="/opt/tools/platanus-v1.2.4/platanus"
exec="/opt/tools/platanus-v1.2.4/platanus_v1.2.4_igonore_quality"
#check if required arguments are there and display usage message
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] || [ -z "$6" ]; then
	printf "\nPlease provide the arguments required for the script.\n\n"
	display_usage
	exit 1	
else
	threads="-t $4"
	memory="-m $5"
	k_size="-k $6"
	kmer_coverage="-c 5"
	output="-o output"
	string_samples="-f"
fi


#generate general command
command_contigs="$exec assemble $output $threads $memory $k_size $kmer_coverage"
command_scaffolds="$exec scaffold $output $threads"
matepairFlag=false



##################################################################
###################Create command#################################
##################################################################
if [ -z "$7" ] && [ -z "$8" ] ; then

	if [[ "$1" != *".f[aq]"* || "$1" != *".fast[aq]"* ]] ; then
 		#contigs 
		assembly_contigs $1 $2 $3
		generate_final_file_contigs 
		#scaffolds
		if [ -z "$9" ]  ; then
		   	scaffolding $1 $2 $3
	            	generate_final_file_scaffolds $8
        	else
			scaffolding $1 $2 $3 $9
            		generate_final_file_scaffolds $8
        	fi

	else	
		#contigs only
		assembly_concatenated_file $1
		printf "\nDone. Note that as only one concatenated fastq file was provided, no scaffolding command was created because it needs the paired end and/or matepair information.\n\n"
        fi

elif [[ -z "$8"  &&  "$7" = "yes" ]] || [[ "$7" = "yes"  &&  "$8" = "no" ]]; then
	if [[ "$1" != *".f[aq]"* || "$1" != *".fast[aq]"* ]] ; then
		assembly_contigs $1 $2 $3
                generate_final_file_contigs
	else 
		assembly_concatenated_file $1
	fi

elif [ "$7" = "yes" ] && [ "$8" = "yes" ] ; then
	printf "\nPlease verify whether you want to perform either only assemble or only scaffolding. If you want to perform only scaffolding please set the option for only assemble to 'no' If you want to perform only assemble don't set the 8th argument at all.\n\n"
	display_usage
	exit 1

elif [ "$8" = "yes" ]; then
	
	if [ -z "$9" ] ; then
		scaffolding $1 $2 $3
		generate_final_file_scaffolds $8
	else 
		if [ -z "${10}" ] || [ "${10}" = "no" ]; then
			scaffolding $1 $2 $3 $9
			generate_final_file_scaffolds $8
		elif [ "${10}" = "yes" ]; then
			scaffolding $1 $2 $3 $8 ${10}
			generate_final_file_scaffolds $8
		else
			printf "\nPlease set a valid value for the tenth argument: [yes|no].\n\n"
		fi
	fi

else 
	printf "\nPlease check if you provided the right parameters [script is not case sensitive].\n\n"
	display_usage
	exit 1
fi

