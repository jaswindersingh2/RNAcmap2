#!/bin/bash

start=`date +%s`

sel_dca=$2

echo $sel_dca

input="$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
input_dir=$(dirname $input)
seq_id=$(basename $(basename $input) | cut -d. -f1)
program_dir=$(dirname $(readlink -f $0))

path_blastn_database=$3
path_infernal_database=$3


mkdir -p $input_dir/${seq_id}_features #&& mkdir -p $input_dir/${seq_id}_outputs
echo ">"$seq_id > $input_dir/${seq_id}_features/$seq_id.fasta
awk -i inplace '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $input 
tail -n1 $input >> $input_dir/${seq_id}_features/$seq_id.fasta
echo "" >> $input_dir/${seq_id}_features/$seq_id.fasta

feature_dir=$input_dir/${seq_id}_features


###### check if aligned homologous sequences file already exists ############
if [ -f $feature_dir/$seq_id.a2m_msa2 ];	then
        echo ""
        echo "==========================================================================="
        echo "    MSA file $feature_dir/$seq_id.a2m_msa2 from Infernal Pipeline already  "
        echo "    exists for query sequence $feature_dir/$seq_id.fasta.                  "
        echo "                                                                           "
        echo "    Delete existing $feature_dir/$seq_id.a2m_msa2 if want to generate new  "
        echo "    alignment file                                                         "
        echo "==========================================================================="
    	echo ""
else

    #################### check if blastn alignment file ready exists ######################
    if [ -f $feature_dir/$seq_id.bla_msa1 ];       then
	    echo ""
	    echo "============================================================================"
	    echo "    MSA-1 file $feature_dir/$seq_id.bla_msa1 from Infernal Pipeline already "
	    echo "    exists for query sequence $feature_dir/$seq_id.fasta.                   "
	    echo "                                                                            "
	    echo "    Delete existing $feature_dir/$seq_id.bla_msa1 if want to generate new   "
	    echo "    alignment file.                                                         "
	    echo "============================================================================"
		echo ""
    else
        echo ""
        echo "==========================================================================================================================="
        echo "      Running BLASTN for first round of homologous sequence search for query sequence $feature_dir/$seq_id.fasta.          "
        echo "      May take 5 mins to few hours depending on sequence length and no. of homologous sequences in database.               "
        echo "==========================================================================================================================="
        echo ""
        blastn -db $path_blastn_database -query $feature_dir/$seq_id.fasta -out $feature_dir/$seq_id.bla_msa1 -evalue 0.001 -num_descriptions 1 -num_threads 8 -line_length 1000 -num_alignments 50000
    fi
			
	if [ $? -eq 0 ]; then
	    echo ""
	    echo "==========================================================="
        echo "      First round of MSA-1 search completed successfully.  "
	    echo "==========================================================="
	    echo ""
	else
        echo ""
        echo "=============================================================================="
        echo "        Error occured while formatting the nt database.                       "
        echo ""
        echo "  Please raise issue at 'https://github.com/jaswindersingh2/RNAcmap2/issues'  "
        echo "=============================================================================="
        echo ""
        exit 1
    fi

	######## reformat the output ################
    echo ""
    echo "=================================================================================================="
    echo "         Converting $feature_dir/$seq_id.bla_msa1 from BLASTN to $feature_dir/$seq_id.sto_msa1.   "
    echo "=================================================================================================="
    echo ""
	$program_dir/utils/parse_blastn_local.pl $feature_dir/$seq_id.bla_msa1 $feature_dir/$seq_id.fasta $feature_dir/$seq_id.aln_msa1
	sed -i 's/\s.*$//' $feature_dir/$seq_id.aln_msa1
	sed -i "s/$seq_id/$seq_id E=0.0/g" $feature_dir/$seq_id.aln_msa1

	$program_dir/utils/seqkit rmdup -n $feature_dir/$seq_id.aln_msa1 > $feature_dir/temp.aln_msa1
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $feature_dir/temp.aln_msa1 | sed '/^$/d' > $feature_dir/$seq_id.aln_msa1

	$program_dir/utils/reformat.pl fas sto $feature_dir/$seq_id.aln_msa1 $feature_dir/$seq_id.sto_msa1


	if [ $? -eq 0 ]; then
	    echo ""
	    echo "=========================================="
        echo "      Converison completed successfully.  "
	    echo "=========================================="
	    echo ""
	else
        echo ""
        echo "==================================================================================================="
        echo "   Error occured while Converting $feature_dir/$seq_id.bla_msa1 to $feature_dir/$seq_id.sto_msa1   "
        echo " "
        echo "   Please raise issue at 'https://github.com/jaswindersingh2/RNAcmap2/issues'                   "
        echo "===================================================================================================="
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
	for i in `awk '{print $2}' $feature_dir/$seq_id.sto_msa1 | head -n5 | tail -n1 | grep -b -o - | sed 's/..$//'`; do sed -i "s/./&-/$i" $feature_dir/$seq_id.db; done

	#########  add reformated ss from last step to .sto file of blastn ##############
	head -n -1 $feature_dir/$seq_id.sto_msa1 > $feature_dir/temp.sto
	echo "#=GC SS_cons                     "`cat $feature_dir/$seq_id.db` > $feature_dir/temp.txt
	cat $feature_dir/temp.sto $feature_dir/temp.txt > $feature_dir/$seq_id.sto_msa1
	echo "//" >> $feature_dir/$seq_id.sto_msa1

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "=================================================================="
        echo "      Consensus Secondary Structure (CSS) generated successfully. "
	    echo "=================================================================="
	    echo ""
	else
        echo ""
        echo "================================================================================"
        echo "             Error occured while generating structure from RNAfold.             "
        echo " "
        echo "  Please raise issue at 'https://github.com/jaswindersingh2/RNAcmap2/issues' "
        echo "================================================================================"
        echo ""
        exit 1
    fi

	######## run infernal ################
    echo ""
    echo "========================================================================================================================"
    echo "      Building Covariance Model from BLASTN alignment (with SS from RNAfold) from $feature_dir/$seq_id.sto_msa1 file.  "
    echo "========================================================================================================================"
    echo ""
	cmbuild --hand -F $feature_dir/$seq_id.cm_msa1 $feature_dir/$seq_id.sto_msa1

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "================================================================================="
        echo "    Covariance Model (CM) built successfully from $feature_dir/$seq_id.sto_msa1. "
	    echo "================================================================================="
	    echo ""
	else
        echo ""
        echo "==============================================================================================="
        echo "     Error occured while building Covariance Model (CM) from cmbuild.           "
        echo " "
        echo "  Please raise issue at 'https://github.com/jaswindersingh2/RNAcmap2/issues'"
        echo "==============================================================================================="
        echo ""
        exit 1
    fi

    echo ""
    echo "======================================================================="
    echo "       Calibrating the Covariance Model $feature_dir/$seq_id.cm_msa1.  "
    echo "======================================================================="
    echo ""
	cmcalibrate $feature_dir/$seq_id.cm_msa1

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "==========================================================="
        echo "    CM calibrated $feature_dir/$seq_id.cm_msa1 successfully.    "
	    echo "==========================================================="
	    echo ""
	else
        echo ""
        echo "================================================================================"
        echo "     Error occured while calibrating $feature_dir/$seq_id.cm_msa1.              "
        echo " "
        echo "  Please raise issue at 'https://github.com/jaswindersingh2/RNAcmap2/issues'    "
        echo "================================================================================"
        echo ""
        exit 1
    fi

    echo ""
    echo "==========================================================================================================================="
    echo "        Second round of homologous sequences search using the calibrated covariance model $feature_dir/$seq_id.cm_msa1.    "
    echo "                 May take 15 mins to few hours for this step.                                                              "
    echo "==========================================================================================================================="
    echo ""
	cmsearch -o $feature_dir/$seq_id.out_msa2 -A $feature_dir/$seq_id.msa_msa2 --cpu 72 --incE 10.0 $feature_dir/$seq_id.cm_msa1 $path_infernal_database

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "==========================================================="
        echo "      Second round of MSA-2 search completed successfully.  "
	    echo "==========================================================="
	    echo ""
	else
        echo ""
        echo "========================================================================================="
        echo "     Error occured during the second round search using CM $feature_dir/$seq_id.cm_msa1. "
        echo " "
        echo "  Please raise issue at 'https://github.com/jaswindersingh2/RNAcmap2/issues'          "
        echo "========================================================================================="
        echo ""
        exit 1
    fi

	######### reformat the alignment without gaps and dashes  ###############
    echo ""
    echo "============================================================================"
    echo "          Reformatting the output alignment $feature_dir/$seq_id.msa_msa2   "
    echo "============================================================================"
    echo ""

	##### check if .msa_msa2 is not empty  #########
	if [[ -s $feature_dir/$seq_id.msa_msa2 ]] 
	  then 
		esl-reformat --replace acgturyswkmbdhvn:................ a2m $feature_dir/$seq_id.msa_msa2 > $feature_dir/temp.a2m_msa2
		sed -i 's/\s.*$//' $feature_dir/temp.a2m_msa2  # remove everything after space
		sed -i "s/$seq_id/$seq_id E=0.0/g" $feature_dir/temp.a2m_msa2   # 
	else 
	  cat $feature_dir/$seq_id.fasta > $feature_dir/temp.a2m_msa2
	  cat $feature_dir/$seq_id.fasta >> $feature_dir/temp.a2m_msa2
	  sed -i '$ s/.$/./' $feature_dir/temp.a2m_msa2
	fi

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "================================================================"
        echo "   Reformatted the $feature_dir/$seq_id.msa_msa2 successfully.  "
	    echo "================================================================"
	    echo ""
	else
        echo ""
        echo "============================================================================================="
        echo "     Error occured during the refomatting the alignment file $feature_dir/$seq_id.msa_msa2.  "
        echo " "
        echo "  Please raise issue at 'https://github.com/jaswindersingh2/RNAcmap2/issues'              "
        echo "============================================================================================="
        echo ""
        exit 1
    fi

	######### remove duplicates sequences from the alignment ###############
    echo ""
    echo "======================================================================="
    echo "          Removing duplicates from the alignment.                      "
    echo "======================================================================="
    echo ""
	$program_dir/utils/seqkit rmdup -s --only-positive-strand $feature_dir/temp.a2m_msa2 > $feature_dir/$seq_id.a2m_msa2
	$program_dir/utils/seqkit rmdup -n $feature_dir/$seq_id.a2m_msa2 > $feature_dir/temp.a2m_msa2
	cp $feature_dir/temp.a2m_msa2 $feature_dir/$seq_id.a2m_msa2

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
        echo "  Please raise issue at 'https://github.com/jaswindersingh2/RNAcmap2/issues'"
        echo "========================================================================================"
        echo ""
        exit 1
    fi

	############# multiline fasta to single line fasta file   #############
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $feature_dir/$seq_id.a2m_msa2 | sed '/^$/d' > $feature_dir/temp.a2m_msa2 
	############# add query sequence at the top of MSA file and consider top 50000 RNAs  #############
    cat $feature_dir/$seq_id.fasta $feature_dir/temp.a2m_msa2 | sed '/^[[:space:]]*$/d' | head -n100000 > $feature_dir/$seq_id.a2m_msa2 
	

fi


############### check Neff-value of $feature_dir/$seq_id.a2m_msa2  alignment  ############
neff=`$program_dir/GREMLIN_CPP/gremlin_cpp -only_neff -alphabet rna -i $feature_dir/$seq_id.a2m_msa2 | grep 'NEFF' | awk '{ print $3}'`
thres=50.0
echo $neff

if (( $(echo "$neff > $thres" |bc -l) )) || [[ ! -s $feature_dir/$seq_id.msa_msa2 ]]; then
	
	echo ""
	echo "=============================================================="
	echo "Neff-value greater than 50 or No-hit found by Infernal search " 
	echo "               Evalute DCA from MSA-2                         "
	echo "=============================================================="
	echo ""

	############### run dca predictors ############
	if [[ $sel_dca = "gremlin" ]]; then

		$program_dir/GREMLIN_CPP/gremlin_cpp -alphabet rna -i $feature_dir/$seq_id.a2m_msa2 -o $feature_dir/$seq_id.dca &> $feature_dir/$seq_id.log_gremlin

	elif [[ $sel_dca = "plmc" ]]; then

		$program_dir/plmc/bin/plmc -c $feature_dir/$seq_id.dca -a -.ACGUNX -le 20 -lh 0.01 -m 50 $feature_dir/$seq_id.a2m_msa2 &> $feature_dir/$seq_id.log_plmc

	elif [[ $sel_dca = "mfdca" ]]; then

		mfdca compute_fn rna $feature_dir/$seq_id.a2m_msa2 --apc --pseudocount 0.5 --verbose &> temp.log

	elif [[ $sel_dca = "plmdca" ]]; then

		plmdca compute_fn rna $feature_dir/$seq_id.a2m_msa2 --max_iterations 500 --num_threads 40 --apc --verbose &> temp.log
		
	fi

else

	echo ""
	echo "=============================="
	echo "  Neff-value less than 50     "
	echo "  Going for of MSA-3 search   "
	echo "=============================="
	echo ""

	$program_dir/utils/reformat.pl a2m sto $feature_dir/$seq_id.a2m_msa2 $feature_dir/$seq_id.sto_msa2
	sed -i 's/#=GF DE/#=GF DE                          E=0.0/g' $feature_dir/$seq_id.sto_msa2  # missing line at the top with E=0.0

	RNAfold $feature_dir/$seq_id.fasta | awk '{print $1}' | tail -n +3 > $feature_dir/$seq_id.db

	################ reformat ss with according to gaps in reference sequence of .sto file from blastn ################
	for i in `awk '{print $2}' $feature_dir/$seq_id.sto_msa2 | head -n5 | tail -n1 | grep -b -o - | sed 's/..$//'`; do sed -i "s/./&-/$i" $feature_dir/$seq_id.db; done

	#########  add reformated ss from last step to .sto file of blastn ##############
	head -n -1 $feature_dir/$seq_id.sto_msa2 > $feature_dir/temp.sto
	echo "#=GC SS_cons                     "`cat $feature_dir/$seq_id.db` > $feature_dir/temp.txt
	cat $feature_dir/temp.sto $feature_dir/temp.txt > $feature_dir/$seq_id.sto_msa2
	echo "//" >> $feature_dir/$seq_id.sto_msa2


	######## run infernal round-2 ################
    echo ""
    echo "==============================================================================================================================="
    echo "      Building Covariance Model from infernal MSA-2 alignment (with SS from RNAfold) from $feature_dir/$seq_id.sto_msa2 file.  "
    echo "==============================================================================================================================="
    echo ""
	cmbuild --hand -F $feature_dir/$seq_id.cm_msa2 $feature_dir/$seq_id.sto_msa2

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "================================================================================="
        echo "    Covariance Model (CM) built successfully from $feature_dir/$seq_id.sto_msa2. "
	    echo "================================================================================="
	    echo ""
	else
        echo ""
        echo "==============================================================================================="
        echo "     Error occured while building Covariance Model (CM) from cmbuild.           "
        echo " "
        echo "  Please raise issue at 'https://github.com/jaswindersingh2/RNAcmap2/issues'"
        echo "==============================================================================================="
        echo ""
        exit 1
    fi

    echo ""
    echo "======================================================================="
    echo "       Calibrating the Covariance Model $feature_dir/$seq_id.cm_msa2.  "
    echo "======================================================================="
    echo ""
	cmcalibrate $feature_dir/$seq_id.cm_msa2

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "==========================================================="
        echo "    CM calibrated $feature_dir/$seq_id.cm_msa2 successfully.    "
	    echo "==========================================================="
	    echo ""
	else
        echo ""
        echo "================================================================================"
        echo "     Error occured while calibrating $feature_dir/$seq_id.cm_msa2.              "
        echo " "
        echo "  Please raise issue at 'https://github.com/jaswindersingh2/RNAcmap2/issues' "
        echo "================================================================================"
        echo ""
        exit 1
    fi

    echo ""
    echo "==========================================================================================================================="
    echo "        Second round of homologous sequences search using the calibrated covariance model $feature_dir/$seq_id.cm_msa2.    "
    echo "                 May take 15 mins to few hours for this step.                                                              "
    echo "==========================================================================================================================="
    echo ""
	cmsearch -o $feature_dir/$seq_id.out_msa3 -A $feature_dir/$seq_id.msa_msa3 --cpu 16 --incE 10.0 $feature_dir/$seq_id.cm_msa2 $path_infernal_database

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "==========================================================="
        echo "      Third round of MSA-3 search completed successfully.  "
	    echo "==========================================================="
	    echo ""
	else
        echo ""
        echo "========================================================================================="
        echo "     Error occured during the second round search using CM $feature_dir/$seq_id.cm_msa2. "
        echo " "
        echo "  Please raise issue at 'https://github.com/jaswindersingh2/RNAcmap2/issues'          "
        echo "========================================================================================="
        echo ""
        exit 1
    fi

	######### reformat the alignment without gaps and dashes  ###############
    echo ""
    echo "============================================================================"
    echo "          Reformatting the output alignment $feature_dir/$seq_id.msa_msa3   "
    echo "============================================================================"
    echo ""

	##### check if .msa_msa3 is not empty  #########
	if [[ -s $feature_dir/$seq_id.msa_msa3 ]] 
	  then 
		esl-reformat --replace acgturyswkmbdhvn:................ a2m $feature_dir/$seq_id.msa_msa3 > $feature_dir/temp.a2m_msa3
	else 
	  cat $feature_dir/$seq_id.fasta > $feature_dir/temp.a2m_msa3
	  cat $feature_dir/$seq_id.fasta >> $feature_dir/temp.a2m_msa3
	  sed -i '$ s/.$/./' $feature_dir/temp.a2m_msa3
	fi

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "================================================================"
        echo "   Reformatted the $feature_dir/$seq_id.msa_msa3 successfully.  "
	    echo "================================================================"
	    echo ""
	else
        echo ""
        echo "============================================================================================="
        echo "     Error occured during the refomatting the alignment file $feature_dir/$seq_id.msa_msa3.  "
        echo " "
        echo "  Please raise issue at 'https://github.com/jaswindersingh2/RNAcmap2/issues'              "
        echo "============================================================================================="
        echo ""
        exit 1
    fi

	######### remove duplicates sequences from the alignment ###############
    echo ""
    echo "======================================================================="
    echo "          Removing duplicates from the alignment.                      "
    echo "======================================================================="
    echo ""
	$program_dir/utils/seqkit rmdup -s --only-positive-strand $feature_dir/temp.a2m_msa3 > $feature_dir/$seq_id.a2m_msa3

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "==============================================="
        echo "   Duplicate sequences removed successfully.   "
	    echo "==============================================="
	    echo ""
	else
        echo ""
        echo "========================================================================================"
        echo "     Error occured during the removel of duplicates from MSA-3.  "
        echo " "
        echo "  Please raise issue at 'https://github.com/jaswindersingh2/RNAcmap2/issues'"
        echo "========================================================================================"
        echo ""
        exit 1
    fi

	############# multiline fasta to single line fasta file   #############
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $feature_dir/$seq_id.a2m_msa3 | sed '/^$/d' > $feature_dir/temp.a2m_msa3 
	############# add query sequence at the top of MSA file and consider top 50000 RNAs  #############
    cat $feature_dir/$seq_id.fasta $feature_dir/temp.a2m_msa3 | sed '/^[[:space:]]*$/d' | head -n100000 > $feature_dir/$seq_id.a2m_msa3 
#	sed -i '/^[[:space:]]*$/d' $feature_dir/$seq_id.a2m_msa3 


	############### run dca predictors ############
	if [[ $sel_dca = "gremlin" ]]; then

		$program_dir/GREMLIN_CPP/gremlin_cpp -alphabet rna -i $feature_dir/$seq_id.a2m_msa3 -o $feature_dir/$seq_id.dca &> $feature_dir/$seq_id.log_gremlin

	elif [[ $sel_dca = "plmc" ]]; then

		$program_dir/plmc/bin/plmc -c $feature_dir/$seq_id.dca -a -.ACGUNX -le 20 -lh 0.01 -m 50 $feature_dir/$seq_id.a2m_msa3 &> $feature_dir/$seq_id.log_plmc

	elif [[ $sel_dca = "mfdca" ]]; then

		mfdca compute_fn rna $feature_dir/$seq_id.a2m_msa3 --apc --pseudocount 0.5 --verbose &> temp.log

	elif [[ $sel_dca = "plmdca" ]]; then

		plmdca compute_fn rna $feature_dir/$seq_id.a2m_msa3 --max_iterations 500 --num_threads 40 --apc --verbose &> temp.log
		
	fi
fi

end=`date +%s`

runtime=$((end-start))

echo -e "\ncomputation time = "$runtime" seconds"

