from Bio import SeqIO
import sys, os
import gzip
from collections import defaultdict
import csv
from sklearn.feature_extraction import DictVectorizer

def loadTerms(inPath, proteins):
    print "Loading terms from", inPath
    with gzip.open(inPath, "rt") as f:
        tsv = csv.reader(f, delimiter='\t')
        for row in tsv:
            protId, goTerm, evCode = row
            protein = proteins[protId]
            if "terms" not in protein:
                protein["terms"] = {}
            protein["terms"][goTerm] = evCode

def loadSequences(inPath, proteins):
    print "Loading sequences from", inPath
    with gzip.open(inPath, "rt") as f:
        sequences = SeqIO.parse(f, 'fasta')
        for seq in sequences:
            proteins[seq.id]["seq"] = seq
            #print seq.id, seq.seq

def buildExamples(proteins):
    print "Building examples"
    for protId in sorted(proteins.keys()):
        print "Processing", protId
        protein = proteins[protId]
        features = {}
        seq = protein["seq"].seq
        for i in range(len(seq)-3):
            feature = seq[i:i+3]
            features[feature] = 1
        print features
        example = DictVectorizer(sparse=False).fit_transform(features)
            

def run(dataPath):
    proteins = defaultdict(lambda: dict())
    loadSequences(os.path.join(options.dataPath, "Swiss_Prot", "Swissprot_sequence.tsv.gz"), proteins)
    loadTerms(os.path.join(options.dataPath, "Swiss_Prot", "Swissprot_evidence.tsv.gz"), proteins)
    print proteins["14310_ARATH"]
    buildExamples(proteins)

if __name__=="__main__":       
    from optparse import OptionParser
    optparser = OptionParser(description="")
    optparser.add_option("-p", "--dataPath", default=os.path.expanduser("~/data/CAFA3"), help="")
    (options, args) = optparser.parse_args()
    
    #proteins = de
    #importProteins(os.path.join(options.dataPath, "Swiss_Prot", "Swissprot_sequence.tsv.gz"))
    run(options.dataPath)