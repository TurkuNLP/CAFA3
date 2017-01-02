from Bio import SeqIO
import sys, os

def importProteins(inPath):
    sequences = SeqIO.parse(open(inPath),'fasta')
    for seq in sequences:
        print seq.id, seq.tostring()

if __name__=="__main__":       
    from optparse import OptionParser
    optparser = OptionParser(description="")
    optparser.add_option("-p", "--dataPath", default=os.path.expanduser("~/data/CAFA3"), help="")
    (options, args) = optparser.parse_args()
    
    importProteins(os.path.join(options.dataPath, "Swiss_Prot"))
    