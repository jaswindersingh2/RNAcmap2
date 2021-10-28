#!/bin/bash




####### environment samples database download ###########
mkdir -p database/env_nt_database

echo ""
echo "==========================================" 
echo " Downloading environment samples database "
echo "==========================================" 
echo ""

curl -o ./database/env_nt_database/env_nt.00.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.00.tar.gz
curl -o ./database/env_nt_database/env_nt.01.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.01.tar.gz
curl -o ./database/env_nt_database/env_nt.02.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.02.tar.gz
curl -o ./database/env_nt_database/env_nt.03.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.03.tar.gz
curl -o ./database/env_nt_database/env_nt.04.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.04.tar.gz
curl -o ./database/env_nt_database/env_nt.05.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.05.tar.gz
curl -o ./database/env_nt_database/env_nt.06.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.06.tar.gz
curl -o ./database/env_nt_database/env_nt.07.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.07.tar.gz
curl -o ./database/env_nt_database/env_nt.08.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.08.tar.gz
curl -o ./database/env_nt_database/env_nt.09.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.09.tar.gz
curl -o ./database/env_nt_database/env_nt.10.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.10.tar.gz
curl -o ./database/env_nt_database/env_nt.11.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.11.tar.gz
curl -o ./database/env_nt_database/env_nt.12.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.12.tar.gz
curl -o ./database/env_nt_database/env_nt.13.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.13.tar.gz
curl -o ./database/env_nt_database/env_nt.14.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.14.tar.gz
curl -o ./database/env_nt_database/env_nt.15.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.15.tar.gz
curl -o ./database/env_nt_database/env_nt.16.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.16.tar.gz
curl -o ./database/env_nt_database/env_nt.17.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.17.tar.gz
curl -o ./database/env_nt_database/env_nt.18.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.18.tar.gz
curl -o ./database/env_nt_database/env_nt.19.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.19.tar.gz
curl -o ./database/env_nt_database/env_nt.20.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.20.tar.gz
curl -o ./database/env_nt_database/env_nt.21.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.21.tar.gz
curl -o ./database/env_nt_database/env_nt.22.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.22.tar.gz
curl -o ./database/env_nt_database/env_nt.23.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.23.tar.gz
curl -o ./database/env_nt_database/env_nt.24.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.24.tar.gz
curl -o ./database/env_nt_database/env_nt.25.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.25.tar.gz
curl -o ./database/env_nt_database/env_nt.26.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.26.tar.gz
curl -o ./database/env_nt_database/env_nt.27.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.27.tar.gz
curl -o ./database/env_nt_database/env_nt.28.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.28.tar.gz
curl -o ./database/env_nt_database/env_nt.29.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.29.tar.gz
curl -o ./database/env_nt_database/env_nt.30.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.30.tar.gz
curl -o ./database/env_nt_database/env_nt.31.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.31.tar.gz
curl -o ./database/env_nt_database/env_nt.32.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.32.tar.gz
curl -o ./database/env_nt_database/env_nt.33.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.33.tar.gz
curl -o ./database/env_nt_database/env_nt.34.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.34.tar.gz
curl -o ./database/env_nt_database/env_nt.35.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/env_nt.35.tar.gz

echo ""
echo "==========================================" 
echo "    Unzip environment samples database    "
echo "==========================================" 
echo ""
cat ./database/env_nt_database/*.tar.gz | tar zxvf - -i -C ./database/env_nt_database/


echo ""
echo "============================================================" 
echo " Get fasta file from formatted environment samples database "
echo "============================================================" 
echo ""
blastdbcmd -db ./database/env_nt_database/env_nt -entry all -outfmt %f -out ./database/env_nt_database/env_nt.fasta






####### transcriptome shotgun assembly database download ###########
mkdir ./database/tsa_nt_database

echo ""
echo "=====================================================" 
echo " Downloading transcriptome shotgun assembly database "
echo "=====================================================" 
echo ""

curl -o ./database/tsa_nt_database/tsa_nt.00.tar.gz "ftp://ftp.ncbi.nlm.nih.gov/blast/db/tsa_nt.00.tar.gz"
curl -o ./database/tsa_nt_database/tsa_nt.01.tar.gz "ftp://ftp.ncbi.nlm.nih.gov/blast/db/tsa_nt.01.tar.gz"
curl -o ./database/tsa_nt_database/tsa_nt.02.tar.gz "ftp://ftp.ncbi.nlm.nih.gov/blast/db/tsa_nt.02.tar.gz"
curl -o ./database/tsa_nt_database/tsa_nt.03.tar.gz "ftp://ftp.ncbi.nlm.nih.gov/blast/db/tsa_nt.03.tar.gz"

echo ""
echo "=====================================================" 
echo "    Unzip transcriptome shotgun assembly database    "
echo "=====================================================" 
echo ""
cat ./database/tsa_nt_database/*.tar.gz | tar zxvf - -i -C ./database/tsa_nt_database/


echo ""
echo "=======================================================================" 
echo " Get fasta file from formatted transcriptome shotgun assembly database "
echo "=======================================================================" 
echo ""
blastdbcmd -db ./database/tsa_nt_database/tsa_nt -entry all -outfmt %f -out ./database/tsa_nt_database/tsa_nt.fasta






#######  Patent Division of GenBank database download ###########
mkdir ./database/pat_nt_database

echo ""
echo "=====================================================" 
echo "   Downloading Patent Division of GenBank database   "
echo "=====================================================" 
echo ""

curl -o ./database/pat_nt_database/patnt.00.tar.gz "ftp://ftp.ncbi.nlm.nih.gov/blast/db/patnt.00.tar.gz"
curl -o ./database/pat_nt_database/patnt.01.tar.gz "ftp://ftp.ncbi.nlm.nih.gov/blast/db/patnt.01.tar.gz"
curl -o ./database/pat_nt_database/patnt.02.tar.gz "ftp://ftp.ncbi.nlm.nih.gov/blast/db/patnt.02.tar.gz"
curl -o ./database/pat_nt_database/patnt.03.tar.gz "ftp://ftp.ncbi.nlm.nih.gov/blast/db/patnt.03.tar.gz"
curl -o ./database/pat_nt_database/patnt.04.tar.gz "ftp://ftp.ncbi.nlm.nih.gov/blast/db/patnt.04.tar.gz"
curl -o ./database/pat_nt_database/patnt.05.tar.gz "ftp://ftp.ncbi.nlm.nih.gov/blast/db/patnt.05.tar.gz"
curl -o ./database/pat_nt_database/patnt.06.tar.gz "ftp://ftp.ncbi.nlm.nih.gov/blast/db/patnt.06.tar.gz"
curl -o ./database/pat_nt_database/patnt.07.tar.gz "ftp://ftp.ncbi.nlm.nih.gov/blast/db/patnt.07.tar.gz"


echo ""
echo "=====================================================" 
echo "      Unzip Patent Division of GenBank database      "
echo "=====================================================" 
echo ""
cat ./database/pat_nt_database/*.tar.gz | tar zxvf - -i -C ./database/pat_nt_database/


echo ""
echo "=======================================================================" 
echo "   Get fasta file from formatted Patent Division of GenBank database   "
echo "=======================================================================" 
echo ""
blastdbcmd -db ./database/pat_nt_database/pat_nt -entry all -outfmt %f -out ./database/pat_nt_database/pat_nt.fasta





#######  Nucleotide database download ###########
mkdir ./database/nt_database

echo ""
echo "=======================================" 
echo "    Downloading Nucleotide (nt) database    "
echo "=======================================" 
echo ""
curl -o ./database/nt_database/nt.gz "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz"


echo ""
echo "=====================================================" 
echo "      Unzip Nucleotide (nt) database      "
echo "=====================================================" 
echo ""
gunzip ./database/nt_database/nt.gz





echo ""
echo "=========================================================" 
echo "    Concatenating four databases into single fasta file  "
echo "=========================================================" 
echo ""
mkdir ./database/combined_database
cat ./database/nt_database/nt ./database/env_nt_database/env_nt.fasta ./database/tsa_nt_database/tsa_nt.fasta ./database/pat_nt_database/pat_nt.fasta > ./database/combined_database/combined.fasta





echo ""
echo "============================================================" 
echo "    Remove duplicate sequences from the combined database   "
echo "============================================================" 
echo ""
seqkit rmdup -s ./database/combined_database/combined.fasta > ./database/combined_database/combined_nr.fasta


