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
        for row in reader:
            if protein == None or protein["id"] != row["id"]:
                assert "predictions" not in protein
                protein["predictions"] = {}
                assert "terms" not in protein
                protein["terms"] = {}
            if row["predicted"] == "1":
                protein["predictions"][row["label"]] = 1
            elif row["gold"] == "1":
                protein["terms"][row["gold"]] = "1"

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
    
    loadPredictions(inPath, proteins)
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
    evaluateFile(options.dataPath, options.output, options.setNames)