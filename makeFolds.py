import os
import loading
from collections import defaultdict
import random
import gzip
import csv

def loadFolds(proteins, inPath, checkSplit=True):
    counts = defaultdict(int)
    print "Loading folds from", inPath
#     with open(inPath, "rt") as f:
#         for line in f:
#             print line
    with gzip.open(inPath, "rt") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row["id"] in proteins:
                counts["row"] += 1
                protein = proteins[row["id"]]
                assert "fold" not in protein
                if checkSplit:
                    assert protein["split"] == row["split"]
                protein["fold"] = int(row["fold"])
            else:
                counts["protein-not-found"] += 1
    print "Loaded folds,", dict(counts)

def saveFolds(proteins, outPath):
    print "Saving folds to", outPath
    outDir = os.path.dirname(outPath)
    if not os.path.exists(outDir):
        print "Making output directory", outDir
        os.makedirs(outDir)
    rows = []
    for key in sorted(proteins.keys()):
        protein = proteins[key]
        assert protein.get("fold") != None
        rows.append({"id":str(key), "fold":str(protein.get("fold")), "split":str(protein.get("split"))})
    with gzip.open(outPath, "wb") as f:
        dw = csv.DictWriter(f, ["id", "fold", "split"], delimiter="\t")
        dw.writeheader()
        dw.writerows(rows)

def makeFolds(proteins, numFolds=10, seed=1):
    print "Generating", str(numFolds) + "-fold division"
    rand = random.Random(seed)
    folds = range(numFolds)
    counts = defaultdict(int)
    for key in sorted(proteins.keys()):
        protein = proteins[key]
        if protein.get("split") == None:
            continue
        protein["fold"] = rand.choice(folds)
        counts[protein["fold"]] += 1
    print "Folds:", dict(counts)
    print "Total:", sum(counts.values())

def run(dataPath, outPath=None):
    print "==========", "Generating Folds", "=========="
    proteins = defaultdict(lambda: dict())
    print "Loading Swissprot proteins"
    loading.loadFASTA(os.path.join(options.dataPath, "Swiss_Prot", "Swissprot_sequence.tsv.gz"), proteins)
    #print "Loading CAFA3 targets"
    #loading.loadFASTA(os.path.join(options.dataPath, "CAFA3_targets", "Target_files", "target.all.fasta"), proteins, True)
    #print "Proteins:", len(proteins)
    #termCounts = loading.loadAnnotations(os.path.join(options.dataPath, "data", "Swissprot_propagated.tsv.gz"), proteins)
    #print "Unique terms:", len(termCounts)
    loading.loadSplit(os.path.join(options.dataPath, "data"), proteins)
    makeFolds(proteins)
    saveFolds(proteins, outPath)

if __name__=="__main__":       
    from optparse import OptionParser
    optparser = OptionParser(description="")
    optparser.add_option("-p", "--dataPath", default=os.path.expanduser("~/data/CAFA3/data"), help="")
    optparser.add_option("-o", "--output", default=None, help="")
    (options, args) = optparser.parse_args()
    
    run(options.dataPath, options.output)