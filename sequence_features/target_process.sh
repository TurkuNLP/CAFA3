#!/bin/bash

main_dir=/home/sukaew/CAFA_test
input_dir=$main_dir/targetFiles
input_seq=$input_dir/sequences

for f in `cat $input_dir/fasta_filelist`
do
    cd $HOME/CAFA3/ncbi-blast-2.5.0+/bin
    export PATH=$PATH:$HOME/CAFA3/ncbi-blast-2.5.0+/bin
    export BLASTDB=$HOME/CAFA3/blastdb
    echo "BLAST62: $f"
    ./blastp -query $f -out $f.xml -outfmt 5 -num_threads 22 -db Swissprot_sequence.tsv -matrix BLOSUM62
    gzip $f.xml
    mv $f.xml.gz $input_dir/blastp62
    echo "BLAST45: $f"
    ./blastp -query $f -out $f.xml -outfmt 5 -evalue 10 -num_threads 22 -db Swissprot_sequence.tsv -matrix BLOSUM45 -max_target_seqs 70000 -gapopen 10 -gapextend 3 -window_size 50
    gzip $f.xml
    mv $f.xml.gz $input_dir/blastp45
    echo "deltaBLAST: $f"
    ./deltablast -query $f -out $f.xml -outfmt 5 -evalue 0.001 -num_threads 20 -db cdd_delta
    gzip $f.xml
    mv $f.xml.gz $input_dir/deltaBlast
    echo "interproscan: $f"
    cd $HOME/CAFA3/interproscan-5.20-59.0
    ./interproscan.sh -i $f -f xml -d $input_seq
    gzip $f.xml
    mv $f.xml.gz $input_dir/interproscan 
done
