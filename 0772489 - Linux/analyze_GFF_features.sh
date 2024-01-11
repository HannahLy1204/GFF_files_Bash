#!/bin/bash

#Checking if input argument is provided
usage="usage: $0 <chromosomeID>. Please provide chromosome ID"
if [ "$#" -lt 1 ]
then
  echo "$usage"
  exit -1
fi

chrID=$1
array=( $(seq 23) ) #Array that contains list of all possible input
array[23]="MT"
array[24]="X"
array[25]="Y"
valid_argument=false

#Checking if input argument is valid
for id in "${array[@]}"; do
	if [ "$id" == "$chrID" ]; then
		valid_argument=true
		break
	fi
done

if ! $valid_argument; then
	echo "$chrID is not a valid chromosomeID (possible value : ${array[@]})"
	exit -1
fi

#Create a result file to store result output
result=result.chr$chrID.txt

#Retrieving data based on the chromosomeID if file is not yet exist in the system
gff3file="Homo_sapiens.GRCh38.110.chromosome.$chrID.gff3"
if [ ! -e $gff3file ]; then
      	wget "https://ftp.ensembl.org/pub/current_gff3/homo_sapiens/$gff3file.gz"
	gunzip "$gff3file.gz"
fi

#Counting features in the gff3 file
echo -e "$gff3file is already downloaded and unzipped. Ready to be analyzed.\n" > $result
echo "Feature count chromosome $chrID" >> $result
echo "-------------------------------" >> $result

#Features of gff3 file lying in column 3 of the file
featurelist=featurelist.chr$chrID.txt
grep "^$chrID" $gff3file|cut -f 3 | sort | uniq -c > $featurelist
cat $featurelist >> $result
echo " " >> $result

#Transcript, gene description (transcript features contains transcriptID and gene_id but no description -> only gene features contains description of the gene)
genelist=transcriptID.geneID.chr$chrID.txt
awk '($3 ~ /transcript$/ || $3 ~ /RNA$/)' $gff3file | awk -F 'transcript:|;Parent=gene:' '{print $2, $3}' | cut -d ';' -f 1| uniq | sort -k 2 > $genelist

description=geneID.description.chr$chrID.txt
awk '($3 ~ /gene$/)' $gff3file |awk -F 'ID=gene:' '{print $2}' |cut -d ';' -f 1,4| awk -F ';description=' '{print $1, $2}'|uniq | sort -k 1 > $description
#Associative array of genes description (gene_id is key for its corresponding description)
declare -A genes
while read -r key desc; do
        genes["$key"]=$desc
done < "$description"

#Check if transcriptID belong to gene_id having the corresponding description
database=transcriptID.database.chr$chrID.txt
while read -r id gene_id; do
    if [ -n "${genes["$gene_id"]}" ]; then
        echo "$id $gene_id ${genes["$gene_id"]}"
    else
        echo "$id $gene_id no description"
    fi
done < "$genelist" | sort -k 1 > "$database"
rm $genelist #remove the old file

#Top 10 transcriptIDs with highest occurences of exons, CDSs, five_prime_UTR, three_prime_UTR
echo "Top 10: chromosome $chrID:" >> $result
echo "--------------------------" >> $result
echo " " >> $result

#Associative array of transcriptID (as key) vs gene_id and description information
declare -A transcriptID
while read -r ids info; do
	transcriptID["$ids"]=$info
done < "$database"

#Display information corresponding to the requirements
my_array=("exon" "CDS" "five_prime_UTR" "three_prime_UTR")
for feature in "${my_array[@]}"
do
	#Counting number of feature for  each transcriptID
        file=$feature.txt
        awk -v var="$feature" '($3 == var)' $gff3file | awk -F'=transcript:' '{print $2}'| cut -d ';' -f 1 |uniq -c | sort -k 1 -n -r > $file
        echo ">>transcriptIDs with the highest number of $feature:"

	#Re-arrange information of $file and write to an output file correspond to the required syntax
	output=$feature.chr$chrID.txt

	#Check if transcriptID available in transcript-gene-description file
	while read -r num trans_id; do
		if [ -n "${transcriptID["$trans_id"]}" ]; then
			echo "Transcript $trans_id >>> #"$feature":$num gene:${transcriptID["$trans_id"]}"
		else
			echo "Transcript $trans_id >>> #"$feature":$num gene:could not find gene_id of this transcript"
		fi
	done < "$file" > "$output"
	rm $file #remove the old file
	head -10 $output
	echo " "
done >> $result 
cat $result
echo "Finish!"
