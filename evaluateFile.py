import os
import loading
from collections import defaultdict
import gzip
import csv
import evaluation

def makeExamples(proteins, limitToSets):
    print "Converting proteins to examples"
    examples = {"labels":[], "ids":[], "cafa_ids":[], "sets":[], "label_names":[], "label_size":{}}
    protIds = sorted(proteins.keys())
    protObjs = [proteins[key] for key in protIds]
    for protein in protObjs:
        if not any(x in limitToSets for x in protein["sets"]):
            continue
        # Build labels
        labels = protein["terms"].keys()
        labels = sorted(labels)
        #if len(labels) == 0:
        #    labels = ["no_annotations"]
        for label in labels:
            if label not in examples["label_size"]:
                examples["label_size"][label] = 0
            examples["label_size"][label] += 1
        examples["labels"].append(labels)
        examples["ids"].append(protein["id"])
        examples["cafa_ids"].append(protein["cafa_ids"])
        examples["sets"].append(protein["sets"])
    print "Converted", len(proteins), "proteins into", len(examples["labels"]), "examples"
    return examples

def loadPredictions(proteins, inPath):
    print "Loading predictions from", inPath
    with gzip.open(inPath, "rt") as f:
        reader = csv.DictReader(f, delimiter='\t')
        protein = None
        counts = {"duplicates":0, "predicted":0, "protein-not-loaded":0}
        skip = False
        for row in reader:
            if protein == None or protein["id"] != row["id"]:
                if row["id"] in proteins:
                    protein = proteins[row["id"]]
                    if "predictions" in protein:
                        counts["duplicates"] += 1
                        skip = True
                    else:
                        skip = False
                        counts["predicted"] += 1
                        protein["predictions"] = {}
                        assert "gold" not in protein
                        protein["gold"] = {}
                else:
                    counts["protein-not-loaded"] += 1
                    skip = True
            if not skip:
                if row["predicted"] == "1":
                    protein["predictions"][row["label"]] = 1
                elif row["gold"] == "1":
                    protein["gold"][row["gold"]] = 1
    print "Predictions loaded:", counts

def evaluateFile(inPath, dataPath, setNames):
    print "==========", "Evaluating", "=========="
    proteins = defaultdict(lambda: dict())
    print "Loading Swissprot proteins"
    loading.loadFASTA(os.path.join(options.dataPath, "Swiss_Prot", "Swissprot_sequence.tsv.gz"), proteins)
    print "Loading CAFA3 targets"
    loading.loadFASTA(os.path.join(options.dataPath, "CAFA3_targets", "Target_files", "target.all.fasta"), proteins, True)
    print "Proteins:", len(proteins)
    #termCounts = loading.loadAnnotations(os.path.join(options.dataPath, "data", "Swissprot_propagated.tsv.gz"), proteins)
    #print "Unique terms:", len(termCounts)
    loading.loadSplit(os.path.join(options.dataPath, "data"), proteins)
    loading.defineSets(proteins, "skip")
    
    loadPredictions(proteins, inPath)
    examples = makeExamples(proteins, setNames)
    loading.vectorizeExamples(examples, None)
    evaluation.evaluate(examples["labels"], examples["predicted"], examples, terms=None, averageOnly=True)

if __name__=="__main__":       
    from optparse import OptionParser
    optparser = OptionParser(description="")
    optparser.add_option("-i", "--input", default=None, help="")
    optparser.add_option("-p", "--dataPath", default=os.path.expanduser("~/data/CAFA3/data"), help="")
    optparser.add_option("-s", "--setNames", default=None, help="")
    (options, args) = optparser.parse_args()
    
    options.setNames = [x.strip() for x in options.setNames.split(",")]
    evaluateFile(options.input, options.dataPath, options.setNames)