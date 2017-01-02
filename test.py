from Bio import SeqIO
import sys, os
import gzip

def importProteins(inPath):
    with gzip.open(inPath, "rt") as f:
        sequences = SeqIO.parse(f, 'fasta')
        for seq in sequences:
            print seq.id, seq.seq

if __name__=="__main__":       
    from optparse import OptionParser
    optparser = OptionParser(description="")
    optparser.add_option("-p", "--dataPath", default=os.path.expanduser("~/data/CAFA3"), help="")
    (options, args) = optparser.parse_args()
    
    importProteins(os.path.join(options.dataPath, "Swiss_Prot", "Swissprot_sequence.tsv.gz"))