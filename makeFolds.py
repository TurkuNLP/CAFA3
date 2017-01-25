import os
import loading
from collections import defaultdict
import random

def makeFolds(proteins):
    rand = random.Random(1)
    for key in sorted(proteins.keys()):
        rand.

def run(dataPath, outDir=None):
    print "==========", "Generating Folds", "=========="
    proteins = defaultdict(lambda: dict())
    print "Loading Swissprot proteins"
    loading.loadFASTA(os.path.join(options.dataPath, "Swiss_Prot", "Swissprot_sequence.tsv.gz"), proteins)
    print "Loading CAFA3 targets"
    loading.loadFASTA(os.path.join(options.dataPath, "CAFA3_targets", "Target_files", "target.all.fasta"), proteins, True)
    print "Proteins:", len(proteins)
    termCounts = loading.loadAnnotations(os.path.join(options.dataPath, "data", "Swissprot_propagated.tsv.gz"), proteins)
    print "Unique terms:", len(termCounts)
    loading.loadSplit(os.path.join(options.dataPath, "data"), proteins)

if __name__=="__main__":       
    from optparse import OptionParser
    optparser = OptionParser(description="")
    optparser.add_option("-p", "--dataPath", default=os.path.expanduser("~/data/CAFA3/data"), help="")
    optparser.add_option("-o", "--output", default=None, help="")
    (options, args) = optparser.parse_args()
    
    run(options.dataPath, options.output)