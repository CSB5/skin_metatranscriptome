#The script takes as the first argument, a string that identifies a cluster e.g. "Cluster 0" in the CD-hit clstr file.
#The second file (INPUT) is the CD-hit clstr file.

#https://stackoverflow.com/questions/8880603/loop-through-an-array-of-strings-in-bash

CLUSTID=$1
INPUT=$2

#split the command in parts that need variable expansion, between double quotes, and parts that do not want this, between single quotes. 

sed -n '/^>'"$CLUSTID"'$/, /^>/{ /^>/!p }' "$INPUT" > temp_clst.txt

head -n 1 temp_clst.txt | cut -f2 -d "|" | awk -F '\\.\\.\\.' '{print $1}' > temp_model_gene_name.txt

#the rest of the members in this cluster, pre-sorted from CD-hit's default output
sed 1d temp_clst.txt | grep -Eo 'cluster_[0-9]*' | uniq > temp_clst_members.txt

awk -v var=`cat temp_model_gene_name.txt` '{print var"\t"$0}' temp_clst_members.txt >> model_to_pangene_cluster_assoc.tsv

rm temp_clst.txt
rm temp_model_gene_name.txt
rm temp_clst_members.txt
