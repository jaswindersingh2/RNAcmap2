#!/bin/bash

start=`date +%s`

input="$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
input_dir=$(dirname $input)
seq_id=$(basename $(basename $input) | cut -d. -f1)
program_dir=$(dirname $(readlink -f $0))

path_blastn_database=/home/jaswinder/Documents/project4/database/nt
path_infernal_database=/home/jaswinder/Documents/project4/database/nt.fasta

#path_blastn_database=$program_dir/nt_database/nt      				# set path to the formatted NCBI's database file without extension 
#path_infernal_database=$program_dir/nt_database/nt					# set path to the NCBI's database database file

mkdir -p $input_dir/${seq_id}_features #&& mkdir -p $input_dir/${seq_id}_outputs
echo ">"$seq_id > $input_dir/${seq_id}_features/$seq_id.fasta
awk -i inplace '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $input 
tail -n1 $input >> $input_dir/${seq_id}_features/$seq_id.fasta
echo "" >> $input_dir/${seq_id}_features/$seq_id.fasta

feature_dir=$input_dir/${seq_id}_features
output_dir=$input_dir/${seq_id}_outputs

#exit 1

if [ ! -f $path_infernal_database ];  then
    echo ""
    echo "========================================================================================"
    echo "            Looks like nt database doesn't exists in the path $path_infernal_database.  "
    echo "            If you want to download the database now, please make sure you have enough  "
    echo "            space in mounted directory and internet connection have enough bandwidth as "
    echo "            file is of size 490 GBs after unzip. It may take forever to download if     "
    echo "                                internet is slow!                                       "
    echo "========================================================================================"
    echo ""

    echo -n "Type 'y' for download or any other key to exit: "    
    read userinput

    if [[ $(echo $userinput | tr '[A-Z]' '[a-z]') == 'y' ]]; then

		echo ""
		echo "=============================================================================================="
		echo "       Downloading NCBI's database form ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz link. "
		echo "                                 May take few hours to download.                              "
		echo "=============================================================================================="
		echo ""
		wget -c "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz" -O "$(dirname "$path_infernal_database")/nt.gz"


		if [[ $? -eq 0 ]]; then 
	        echo ""
	        echo "======================================================================="
	        echo "            nt database is completed successfully.                     "
	        echo "======================================================================="
	        echo ""
		else
	        echo ""
	        echo "======================================================================="
	        echo "            Error! Unable to download database sucessfully.            "
	        echo "            Check wget command or internet connection.            "
	        echo "======================================================================="
	        echo ""
	        exit 1        
		fi

		echo ""
		echo "======================================================================"
		echo "            Unziping the downloaded nt database.                      "
		echo "       May take few hours as size of unzipped file is around 270 GBs. "
		echo "======================================================================"
		echo ""
		
	############ unzip the nt data base file ############
		gunzip "$(dirname "$path_infernal_database")/nt.gz"

		if [[ $? -eq 0 ]]; then 
	        echo ""
	        echo "======================================================================="
	        echo "            nt database unzip completed successfully.                  "
	        echo "======================================================================="
	        echo ""
		else
	        echo ""
	        echo "======================================================================="
	        echo "            Error! unable to unzip database sucessfully.               "
	        echo "            Please check if gunzip program exists!                     "
	        echo "======================================================================="
	        echo ""
	        exit 1        
		fi

    else
		echo ""
		echo "==========================================================="
		echo "      Exiting the program because nt database is missing! "
		echo "==========================================================="
		echo ""
        exit 1
    fi

fi


###### check if aligned homologous sequences file already exists ############
if [ -f $feature_dir/$seq_id.a2m ];	then
        echo ""
        echo "======================================================================"
        echo "    MSA file $feature_dir/$seq_id.a2m from Infernal Pipeline already  "
        echo "    exists for query sequence $feature_dir/$seq_id.fasta.             "
        echo "                                                                      "
        echo "    Delete existing $feature_dir/$seq_id.a2m if want to generate new  "
        echo "    alignment file                                                    "
        echo "======================================================================"
    	echo ""
else

   #### check if formatted nt database exists or not ##### 
    if [[ ! -f "$path_blastn_database.nal" ]]; then
        echo ""
        echo "====================================================================="
        echo "    Nucleotide database file $path_blastn_database/nt need to formated      "
        echo "    formated to use with 'makeblastdb' program in BLAST-N program.   "  
        echo ""          
		echo "    Formatting may take several hours as size of file is around 490 GBs. "
        echo "====================================================================="
        echo ""
        makeblastdb -in "$path_infernal_database" -dbtype nucl -out "$path_blastn_database"
        
        if [[ $? -eq 0 ]]; then
                echo ""
                echo "======================================================="
                echo "          nt database formatted successfully.          "
                echo "======================================================="
                echo ""
        else
                echo ""
                echo "=================================================================="
                echo "        Error occured while formatting the nt database.           "
                echo ""
                echo "  Please raise issue at 'https://github.com/jaswindersingh2/SPOT-RNA-2D/issues'"
                echo "=================================================================="
                echo ""
                exit 1
        fi                      
    fi


    #################### check if blastn alignment file ready exists ######################
    if [ -f $feature_dir/$seq_id.bla ];       then
	    echo ""
	    echo "======================================================================="
	    echo "    MSA-1 file $feature_dir/$seq_id.bla from Infernal Pipeline already "
	    echo "    exists for query sequence $feature_dir/$seq_id.fasta.              "
	    echo "                                                                       "
	    echo "    Delete existing $feature_dir/$seq_id.a2m if want to generate new   "
	    echo "    alignment file.                                                    "
	    echo "======================================================================="
		echo ""
    else
        echo ""
        echo "==========================================================================================================================="
        echo "      Running BLASTN for first round of homologous sequence search for query sequence $feature_dir/$seq_id.fasta.          "
        echo "      May take 5 mins to few hours depending on sequence length and no. of homologous sequences in database.               "
        echo "==========================================================================================================================="
        echo ""
        blastn -db $path_blastn_database -query $feature_dir/$seq_id.fasta -out $feature_dir/$seq_id.bla -evalue 0.001 -num_descriptions 1 -num_threads 8 -line_length 1000 -num_alignments 50000
    fi
			
	if [ $? -eq 0 ]; then
	    echo ""
	    echo "==========================================================="
        echo "      First round of MSA-1 search completed successfully.  "
	    echo "==========================================================="
	    echo ""
	else
        echo ""
        echo "=================================================================="
        echo "        Error occured while formatting the nt database.           "
        echo ""
        echo "  Please raise issue at 'https://github.com/jaswindersingh2/SPOT-RNA-2D/issues'"
        echo "=================================================================="
        echo ""
        exit 1
    fi

	######## reformat the output ################
    echo ""
    echo "========================================================================================"
    echo "         Converting $feature_dir/$seq_id.bla from BLASTN to $feature_dir/$seq_id.sto.   "
    echo "========================================================================================"
    echo ""
	$program_dir/utils/parse_blastn_local.pl $feature_dir/$seq_id.bla $feature_dir/$seq_id.fasta $feature_dir/$seq_id.aln
	$program_dir/utils/reformat.pl fas sto $feature_dir/$seq_id.aln $feature_dir/$seq_id.sto


	if [ $? -eq 0 ]; then
	    echo ""
	    echo "=========================================="
        echo "      Converison completed successfully.  "
	    echo "=========================================="
	    echo ""
	else
        echo ""
        echo "============================================================================================="
        echo "   Error occured while Converting $feature_dir/$seq_id.bla to $feature_dir/$seq_id.sto       "
        echo " "
        echo "   Please raise issue at 'https://github.com/jaswindersingh2/SPOT-RNA-2D/issues'             "
        echo "============================================================================================="
        echo ""
        exit 1
    fi

	######## predict secondary structure from RNAfold ################
    echo ""
    echo "==============================================================================================================================="
    echo "       Predicting Consensus Secondary Structure (CSS) of query sequence $feature_dir/$seq_id.fasta using RNAfold predictor.   "
    echo "==============================================================================================================================="
    echo ""

	RNAfold $feature_dir/$seq_id.fasta | awk '{print $1}' | tail -n +3 > $feature_dir/$seq_id.db

	################ reformat ss with according to gaps in reference sequence of .sto file from blastn ################
	for i in `awk '{print $2}' $feature_dir/$seq_id.sto | head -n5 | tail -n1 | grep -b -o - | sed 's/..$//'`; do sed -i "s/./&-/$i" $feature_dir/$seq_id.db; done

	#########  add reformated ss from last step to .sto file of blastn ##############
	head -n -1 $feature_dir/$seq_id.sto > $feature_dir/temp.sto
	echo "#=GC SS_cons                     "`cat $feature_dir/$seq_id.db` > $feature_dir/temp.txt
	cat $feature_dir/temp.sto $feature_dir/temp.txt > $feature_dir/$seq_id.sto
	echo "//" >> $feature_dir/$seq_id.sto

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "=================================================================="
        echo "      Consensus Secondary Structure (CSS) generated successfully. "
	    echo "=================================================================="
	    echo ""
	else
        echo ""
        echo "=============================================================================="
        echo "             Error occured while generating structure from RNAfold.          "
        echo " "
        echo "  Please raise issue at 'https://github.com/jaswindersingh2/SPOT-RNA-2D/issues'"
        echo "=============================================================================="
        echo ""
        exit 1
    fi

	######## run infernal ################
    echo ""
    echo "=============================================================================================================="
    echo "      Building Covariance Model from BLASTN alignment (with SS from SPOT-RNA) from $feature_dir/$seq_id.sto file.         "
    echo "=============================================================================================================="
    echo ""
	cmbuild --hand -F $feature_dir/$seq_id.cm $feature_dir/$seq_id.sto

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "============================================================================"
        echo "    Covariance Model (CM) built successfully from $feature_dir/$seq_id.sto. "
	    echo "============================================================================"
	    echo ""
	else
        echo ""
        echo "==============================================================================================="
        echo "     Error occured while building Covariance Model (CM) from cmbuild.           "
        echo " "
        echo "  Please raise issue at 'https://github.com/jaswindersingh2/SPOT-RNA-2D/issues'"
        echo "==============================================================================================="
        echo ""
        exit 1
    fi

    echo ""
    echo "===================================================================="
    echo "       Calibrating the Covariance Model $feature_dir/$seq_id.cm.    "
    echo "===================================================================="
    echo ""
	cmcalibrate $feature_dir/$seq_id.cm

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "==========================================================="
        echo "    CM calibrated $feature_dir/$seq_id.cm successfully.    "
	    echo "==========================================================="
	    echo ""
	else
        echo ""
        echo "==============================================================="
        echo "     Error occured while calibrating $feature_dir/$seq_id.cm.  "
        echo " "
        echo "  Please raise issue at 'https://github.com/jaswindersingh2/SPOT-RNA-2D/issues'"
        echo "==============================================================="
        echo ""
        exit 1
    fi

    echo ""
    echo "======================================================================================================================"
    echo "        Second round of homologous sequences search using the calibrated covariance model $feature_dir/$seq_id.cm.    "
    echo "                 May take 15 mins to few hours for this step.                                                         "
    echo "======================================================================================================================"
    echo ""
	cmsearch -o $feature_dir/$seq_id.out -A $feature_dir/$seq_id.msa --cpu 24 --incE 10.0 $feature_dir/$seq_id.cm $path_infernal_database

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "==========================================================="
        echo "      Second round of MSA-2 search completed successfully.  "
	    echo "==========================================================="
	    echo ""
	else
        echo ""
        echo "===================================================================================="
        echo "     Error occured during the second round search using CM $feature_dir/$seq_id.cm. "
        echo " "
        echo "  Please raise issue at 'https://github.com/jaswindersingh2/SPOT-RNA-2D/issues'"
        echo "===================================================================================="
        echo ""
        exit 1
    fi

	######### reformat the alignment without gaps and dashes  ###############
    echo ""
    echo "======================================================================="
    echo "          Reformatting the output alignment $feature_dir/$seq_id.msa   "
    echo "          for PSSM and DCA features by removing the gaps and dashes.   "
    echo "======================================================================="
    echo ""

	##### check if .msa	is not empty  #########
	if [[ -s $feature_dir/$seq_id.msa ]] 
	  then 
		esl-reformat --replace acgturyswkmbdhvn:................ a2m $feature_dir/$seq_id.msa > $feature_dir/temp.a2m
	else 
	  cat $feature_dir/$seq_id.fasta > $feature_dir/temp.a2m
	  cat $feature_dir/$seq_id.fasta >> $feature_dir/temp.a2m
	  sed -i '$ s/.$/./' $feature_dir/temp.a2m
	fi

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "==========================================================="
        echo "   Reformatted the $feature_dir/$seq_id.msa successfully.  "
	    echo "==========================================================="
	    echo ""
	else
        echo ""
        echo "========================================================================================"
        echo "     Error occured during the refomatting the alignment file $feature_dir/$seq_id.msa.  "
        echo " "
        echo "  Please raise issue at 'https://github.com/jaswindersingh2/SPOT-RNA-2D/issues'"
        echo "========================================================================================"
        echo ""
        exit 1
    fi

	######### remove duplicates sequences from the alignment ###############
    echo ""
    echo "======================================================================="
    echo "          Removing duplicates from the alignment.                      "
    echo "======================================================================="
    echo ""
	$program_dir/utils/seqkit rmdup -s $feature_dir/temp.a2m > $feature_dir/$seq_id.a2m

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "==============================================="
        echo "   Duplicate sequences removed successfully.   "
	    echo "==============================================="
	    echo ""
	else
        echo ""
        echo "========================================================================================"
        echo "     Error occured during the removel of duplicates from MSA-2.  "
        echo " "
        echo "  Please raise issue at 'https://github.com/jaswindersingh2/SPOT-RNA-2D/issues'"
        echo "========================================================================================"
        echo ""
        exit 1
    fi

	############# multiline fasta to single line fasta file   #############
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $feature_dir/$seq_id.a2m | sed '/^$/d' > $feature_dir/temp.a2m 
	############# add query sequence at the top of MSA file  #############
    cat $feature_dir/$seq_id.fasta $feature_dir/temp.a2m > $feature_dir/$seq_id.a2m 

fi

############# check if pssm file already exists otherwise generate from alignment file #############
if [ -f $feature_dir/$seq_id.pssm ];	then
        echo ""
        echo "=============================================================================================================================================="
        echo "    PSSM feature file $feature_dir/$seq_id.pssm already exists for query sequence $feature_dir/$seq_id.fasta.  "
        echo "=============================================================================================================================================="
    	echo ""
else
	echo ""
	echo "======================================================================================"
	echo "          Extracting PSSM features from the alignment $feature_dir/$seq_id.a2m.       "
	echo "======================================================================================"
	echo ""
	$program_dir/utils/getpssm.pl $feature_dir/$seq_id.fasta $feature_dir/$seq_id.a2m $feature_dir/$seq_id.pssm

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "==============================================================="
        echo "   PSSM extracted successfully from $feature_dir/$seq_id.a2m.  "
	    echo "==============================================================="
	    echo ""
	else
        echo ""
        echo "========================================================================="
        echo "     Error occured while extracting PSSM from $feature_dir/$seq_id.a2m.  "
        echo " "
        echo "  Please raise issue at 'https://github.com/jaswindersingh2/SPOT-RNA-2D/issues'"
        echo "========================================================================="
        echo ""
        exit 1
    fi
fi

############# check if dca file already exists otherwise generate from alignment file #############
if [ -f $feature_dir/$seq_id.dca ];	then
        echo ""
        echo "==============================================================="
        echo "    GRELMLIN feature file $feature_dir/$seq_id.dca already     "
        echo "    exists for query sequence $feature_dir/$seq_id.fasta.      "
        echo " "
        echo "    Delete the existing file if want to generate new dca file. "
        echo "==============================================================="
    	echo ""
else
	echo ""
	echo "============================================================================"
	echo "          Running PLMC for DCA features.                                 "
	echo "============================================================================"
	echo ""
	$program_dir/plmc/bin/plmc -c $feature_dir/$seq_id.dca -a -.ACGUNX -le 20 -lh 0.01 -m 50 $feature_dir/$seq_id.a2m &> $feature_dir/$seq_id.log_plmc
	if [ $? -eq 0 ]; then
		echo ""
		echo "===================================================="
		echo "   DCA features successfully obtained from PLMC. "
		echo "===================================================="
		echo ""
	else
		echo ""
		echo "============================================================================="
		echo "                Error occured while running PLMC.  "
		echo " "
		echo "  Please raise issue at 'https://github.com/jaswindersingh2/SPOT-RNA-2D/issues'"
		echo "============================================================================="
		echo ""
		exit 1
	fi
fi

cp $input_dir/${seq_id}_features/$seq_id.fasta $input_dir/${seq_id}_features/$seq_id

end=`date +%s`

runtime=$((end-start))

echo -e "\ncomputation time = "$runtime" seconds"

