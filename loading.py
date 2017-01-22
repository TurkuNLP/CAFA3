import gzip
import csv
import os
from collections import defaultdict
from sklearn.feature_extraction import DictVectorizer
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.feature_selection.variance_threshold import VarianceThreshold

def openAny(inPath, mode):
    return gzip.open(inPath, mode) if inPath.endswith(".gz") else open(inPath, mode)

def loadAnnotations(inPath, proteins):
    print "Loading annotations from", inPath
    counts = defaultdict(int)
    with gzip.open(inPath, "rt") as f:
        tsv = csv.reader(f, delimiter='\t')
        for row in tsv:
            protId, goTerm, evCode = row
            protein = proteins[protId]
            if "terms" not in protein:
                protein["terms"] = {}
            protein["terms"][goTerm] = evCode
            counts[goTerm] += 1
    return counts

def loadGOTerms(inPath):
    print "Loading GO terms from", inPath
    terms = {}
    with open(inPath, "rt") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            #print row
            assert row["id"] not in terms
            terms[row["id"]] = row
    return terms

def addProtein(proteins, protId, cafaId, sequence, filename, replaceSeq=False, verbose=False, counts=None):
    assert len(sequence) > 0
    if protId not in proteins:
        proteins[protId]["seq"] = sequence
        proteins[protId]["id"] = protId
        proteins[protId]["cafa_ids"] = [cafaId] if cafaId else []
        proteins[protId]["file"] = [filename]
        proteins[protId]["terms"] = {}
        counts["unique"] += 1
    else:
        counts["multiple"] += 1
        if proteins[protId]["seq"] != sequence:
            if verbose:
                print "WARNING, sequence mismatch for", protId, cafaId, ":", (proteins[protId], (protId, cafaId, sequence, filename))
            counts["mismatch"] += 1
            counts["mismatches"].append(protId)
            if replaceSeq:
                proteins[protId]["seq"] = sequence
        assert proteins[protId]["id"] == protId, (proteins[protId], (protId, cafaId, sequence, filename))
        proteins[protId]["file"] += [filename]
        if cafaId != None:
            if len(proteins[protId]["cafa_ids"]) > 0:
                if verbose:
                    print "WARNING, duplicate CAFA target", protId, cafaId, ":", (proteins[protId], (protId, cafaId, sequence, filename))
                counts["duplicate_id"] += 1
                counts["duplicate_ids"].append(protId)
                proteins[protId]["cafa_ids"].append(cafaId)
            else:
                proteins[protId]["cafa_ids"] += [cafaId]
         
def loadFASTA(inPath, proteins, cafaHeader=False):
    print "Loading sequences from", inPath
    filename = os.path.basename(inPath)
    counts = defaultdict(int)
    counts["mismatches"] = []
    counts["duplicate_ids"] = []
    with openAny(inPath, "rt") as f:
        protId = None
        cafaId = None
        sequence = ""
        for line in f:
            if line.startswith(">"):
                # Add already read protein
                if protId != None:
                    addProtein(proteins, protId, cafaId, sequence, filename, counts=counts)
                protId = None
                cafaId = None
                sequence = ""
                # Begin new protein
                protId = line[1:].strip()
                if cafaHeader:
                    cafaId, protId = protId.split()
            else:
                sequence += line.strip()
        if protId != None:
            addProtein(proteins, protId, cafaId, sequence, filename, counts=counts)               
            #print seq.id, seq.seq
    print dict(counts)

def loadSplit(inPath, proteins):
    for dataset in ("train", "devel", "test"):
        filePath = os.path.join(inPath, dataset + ".txt.gz")
        assert os.path.exists(filePath), filePath
        with gzip.open(filePath, "rt") as f:
            for line in f:
                protId = line.strip()
                assert protId in proteins
                proteins[protId]["split"] = dataset

def defineSets(proteins, cafaTargets):
    counts = defaultdict(int)
    for protein in proteins.values():
        cafaSet = ["cafa"] if len(protein["cafa_ids"]) > 0 else []
        splitSet = [protein["split"]] if protein.get("split") != None else []
        if len(cafaSet) > 0:
            if cafaTargets == "overlap":
                protein["sets"] = cafaSet + splitSet
            elif cafaTargets == "separate":
                if "train" in splitSet:
                    protein["sets"] = cafaSet
                else:
                    protein["sets"] = cafaSet + splitSet
            elif cafaTargets == "external":
                protein["sets"] = splitSet if len(splitSet) > 0 else cafaSet
            else:
                raise Exception("CAFA targets were loaded with mode '" + cafaTargets + "'")
        else:
            protein["sets"] = splitSet
        assert len(protein["sets"]) > 0
        category = ",".join(cafaSet + splitSet) + "=>" + ",".join(protein["sets"])
        counts[category] += 1
    print "Defined sets:", dict(counts)

def saveFeatureNames(names, outPath):
    print "Saving feature names to", outPath
    with open(outPath, "wt") as f:
        f.write("index\tname\n")
        for i in range(len(names)):
            f.write(str(i) + "\t" + names[i] + "\n")
    
def vectorizeExamples(examples, featureGroups):
    mlb = MultiLabelBinarizer()
    examples["labels"] = mlb.fit_transform(examples["labels"])
    examples["label_names"] = mlb.classes_
    dv = DictVectorizer(sparse=True)
    examples["features"] = dv.fit_transform(examples["features"])
    examples["feature_names"] = dv.feature_names_
    if featureGroups != None and "select" in featureGroups:
        threshold = .1
        print "Selecting features", examples["features"].shape[1]
        examples["features"] = VarianceThreshold(threshold * (1 - threshold)).fit_transform(examples["features"])
        print "Selected features", examples["features"].shape[1]
        #examples["features"] = SelectKBest(chi2, k=1000).fit_transform(examples["features"], examples["labels"])
    print "Vectorized", len(examples["labels"]), "examples with", len(examples["feature_names"]), "unique features and", len(examples["label_names"]), "unique labels"
