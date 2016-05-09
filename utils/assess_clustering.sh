# This script assess clustering made by pangmer under
# given parameters
#
# It extracts the sequence representative of each cluster
# and all sequences that belongs to that cluster
# Then aligns each of these sequence to the representative
# one using emboss supermatcher. In that way we can obtain
# the similarity for each sequence of the cluster to the
# representative one.
# Finally supermatcher output files are parsed to get
# a summarized file of similarities for every cluster of
# the pangenome
#
#
 
pangenome=$1
mapping=$2
fasta=$3
minlength=$4

# Count number of clusters in the pangenome
ngroups=$(grep '>' $pangenome | wc -l)

# Extract seqs for each cluster
for i in $(seq 1 $ngroups)
  do
    echo Extract representative sequence of the cluster
    Script-fasta-regex-picker -f $pangenome "Group:$i\b" > Group$i.fasta
    echo Count sequence length
    seqlength=$(grep -v ">" Group$i.fasta | wc -m)
    echo "seqlength = $seqlength"
    # Only follow if seqlength is >= to minlength
    if [ $seqlength -ge $minlength ]
      then
         echo Extract sequences forming that cluster
         pang-util-group_seqs -f $fasta -m $mapping -g $i -o Group$i.seqs
         echo Make alignments using supermatcher
          supermatcher -asequence Group$i.fasta -bsequence Group$i.seqs \
          -gapopen 10 -gapextend 0.5 -outfile Group$i.align
	echo "Group$i" >> cluster.similarities
        grep Similarity Group$i.align >> cluster.similarities
        rm Group$i*  
     fi
  done


