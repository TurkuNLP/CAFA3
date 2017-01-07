import sys, os
import gzip
from collections import defaultdict
import csv
from sklearn.feature_extraction import DictVectorizer
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.grid_search import GridSearchCV, ParameterGrid
from sklearn.ensemble.forest import ExtraTreesClassifier, RandomForestClassifier
from sklearn.metrics import classification_report, f1_score
from featureBuilders import *
from utils import Stream
import operator
import time
from sklearn.cross_validation import train_test_split
from sklearn.metrics.ranking import roc_auc_score
from sklearn.metrics.classification import precision_recall_fscore_support
import shutil
import pickle

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

def loadSequences(inPath, proteins):
    print "Loading sequences from", inPath
    with gzip.open(inPath, "rt") as f:
        header = None
        for line in f:
            if header == None:
                assert line.startswith(">")
                header = line[1:].strip()
            else:
                proteins[header]["seq"] = line.strip()
                proteins[header]["id"] = header
                header = None
            #print seq.id, seq.seq

def loadSplit(inPath, proteins):
    for dataset in ("train", "devel", "test"):
        filePath = os.path.join(inPath, dataset + ".txt.gz")
        assert os.path.exists(filePath), filePath
        with gzip.open(filePath, "rt") as f:
            for line in f:
                protId = line.strip()
                assert protId in proteins
                proteins[protId]["set"] = dataset

def saveFeatureNames(names, outPath):
    print "Saving feature names to", outPath
    with open(outPath, "wt") as f:
        f.write("index\tname\n")
        for i in range(len(names)):
            f.write(str(i) + "\t" + names[i] + "\n")
        
# def splitProteins(proteins):
#     datasets = {"devel":[], "train":[], "test":[]}
#     for protId in sorted(proteins.keys()):
#         datasets[proteins[protId]["set"]].append(proteins[protId])
#     print "Divided sets", [(x, len(datasets[x])) for x in sorted(datasets.keys())]
#     return datasets

def buildExamples(proteins, dataPath, limit=None, limitTerms=None, featureGroups=None):
    print "Building examples"
    examples = {"labels":[], "features":None, "ids":[], "sets":[], "label_names":[], "label_size":{}}
    protIds = sorted(proteins.keys())
    if limit:
        protIds = protIds[0:limit]
    protObjs = [proteins[key] for key in protIds]
    for protein in protObjs:
        # Initialize features
        protein["features"] = {"dummy":1}
        # Build labels
        labels = protein["terms"].keys()
        if limitTerms:
            labels = [x for x in labels if x in limitTerms]
        labels = sorted(labels)
        if len(labels) == 0:
            labels = ["no_annotations"]
        for label in labels:
            if label not in examples["label_size"]:
                examples["label_size"][label] = 0
            examples["label_size"][label] += 1
        examples["labels"].append(labels)
        examples["ids"].append(protein["id"])
        examples["sets"].append(protein["set"])
    # Build features
    if featureGroups == None or "similar" in featureGroups:
        builder = UniprotFeatureBuilder(os.path.join(dataPath, "Uniprot", "similar.txt"))
        builder.build(protObjs)
    if featureGroups == None or "blast" in featureGroups:
        builder = BlastFeatureBuilder(os.path.join(dataPath, "blastp_result_features"))
        builder.build(protObjs)
    builder = None
    # Prepare the examples
    mlb = MultiLabelBinarizer()
    dv = DictVectorizer(sparse=True)
    examples["features"] = dv.fit_transform([x["features"] for x in protObjs])
    examples["feature_names"] = dv.feature_names_
    examples["labels"] = mlb.fit_transform(examples["labels"])
    examples["label_names"] = mlb.classes_
    print "Built", len(examples["labels"]), "examples with", len(examples["feature_names"]), "unique features"
    return examples
    #return mlb.fit_transform(examples["labels"]), dv.fit_transform(examples["features"])

def getTopTerms(counts, num=1000):
    return sorted(counts.items(), key=operator.itemgetter(1), reverse=True)[0:num]

# def getResults(examples, scores, terms=None):
#     assert len(scores) == len(examples["label_names"])
#     results = []
#     for i in range(len(examples["label_names"])):
#         label = examples["label_names"][i]
#         result = {"score":scores[i], "id":label, "label_size":examples["label_size"][label]}
#         result["ns"] = None
#         result["name"] = None
#         if terms != None and label in terms:
#             term = terms[label]
#             result["ns"] = term["ns"]
#             result["name"] = term["name"]
#         results.append(result)
#     return results

def metricsToString(result):
    s = "a/f|p/r|s = " 
    s += "%.2f" % result["auc"] + "/" + "%.2f" % result["fscore"]
    s += "|" + "%.2f" % result["precision"] + "/" + "%.2f" % result["recall"] + "|" + "%.2f" % result["support"]
    return s

def printResults(results, maxNumber=None):
    count = 0
    results = sorted(results, key=lambda x: x["auc"], reverse=True)
    for result in sorted(results, reverse=True):
        if result == "average":
            continue
        print metricsToString(result), [result.get("id"), result.get("ns"), result.get("label_size"), result.get("name")]
        count += 1
        if count > maxNumber:
            break

def saveResults(results, outPath):
    print "Writing results to", outPath
    with open(outPath, "wt") as f:
        dw = csv.DictWriter(f, ["auc", "fscore", "precision", "recall", "support", "id", "label_size", "ns", "name"], delimiter='\t')
        dw.writeheader()
        dw.writerow(results["average"])
        results = [x for x in results.values() if x["id"] != "average"]
        dw.writerows(sorted(results, key=lambda x: x["auc"], reverse=True))

def evaluate(labels, predicted, label_names, label_size=None, terms=None):
    results = {}
    # Get the average result
    auc = roc_auc_score(labels, predicted, average="micro")
    p, r, f, s = precision_recall_fscore_support(labels, predicted, average="micro")
    results["average"] = {"id":"average", "ns":None, "name":None, "auc":auc, "precision":p, "recall":r, "fscore":f, "support":s}
    # Get results per label
    aucs = roc_auc_score(labels, predicted, average=None)
    scores = precision_recall_fscore_support(labels, predicted, average=None, labels=label_names)
    assert len(aucs) == len(scores) == len(label_names), (len(aucs), len(scores), len(label_names))
    for auc, score, label_name in zip(aucs, scores, label_names):
        assert label_name not in results
        result = {"id":label_name, "ns":None, "name":None, "auc":auc, 
                  "precision":score[0], "recall":score[1], "fscore":score[2], "support":score[3]}
        results[label_name] = result
        if label_size != None:
            result["label_size"] = label_size[label_name]
        if terms != None and label_name in terms:
            term = terms[label_name]
            result["ns"] = term["ns"]
            result["name"] = term["name"]   
    return results

def optimize(examples, verbose=3, n_jobs = -1, scoring = "f1_micro", cvJobs=1, terms=None):
    grid = ParameterGrid({"n_estimators":[10], "n_jobs":[n_jobs], "verbose":[verbose]}) #{"n_estimators":[1,2,10,50,100]}
    #XTrainAndDevel, XTest, yTrainAndDevel, yTest = train_test_split(X, y, test_size=0.2, random_state=0)
    #XTrain, XDevel, yTrain, yDevel = train_test_split(XTrainAndDevel, yTrainAndDevel, test_size=0.2, random_state=0)
    sets = examples["sets"]
    trainIndices = [i for i in range(len(sets)) if sets[i] == "train"]
    develIndices = [i for i in range(len(sets)) if sets[i] == "devel"]
    trainFeatures = examples["features"][trainIndices]
    develFeatures = examples["features"][develIndices]
    trainLabels = examples["labels"][trainIndices]
    develLabels = examples["labels"][develIndices]
    print "Optimizing, train / devel = ", trainFeatures.shape[0], "/", develFeatures.shape[0]
    best = None
    print "Parameter grid search", time.strftime('%X %x %Z')
    for args in grid:
        print "Learning with args", args
        cls = RandomForestClassifier(**args)
        cls.fit(trainFeatures, trainLabels)
        predicted = cls.predict(develFeatures)
        #score = roc_auc_score(develLabels, predicted, average="micro")
        #scores = roc_auc_score(develLabels, predicted, average=None)
        #print "Average =", score
        #results = getResults(examples, scores, terms)
        results = evaluate(develLabels, predicted, examples["label_names"], examples["label_size"], terms)
        print "Average:", metricsToString(results["average"])
        printResults(results, 20)
        if best == None or results["average"]["auc"] > best["results"]["average"]["auc"]:
            best = {"results":results, "args":args}
        print time.strftime('%X %x %Z')
    return best
    #clf = GridSearchCV(RandomForestClassifier(), args, verbose=verbose, n_jobs=cvJobs, scoring=scoring)
    #clf.fit(X, y)
    #print "Best params", (clf.best_params_, clf.best_score_)

# def learn(train, devel, test, limit=None, limitTerms=None, featureGroups=None):
#     print time.strftime('%X %x %Z')
#     print "Building devel examples"
#     develLabels, develFeatures = buildExamples(devel, limit, limitTerms, featureGroups)
#     print "Building train examples"
#     trainLabels, trainFeatures = buildExamples(train, limit, limitTerms, featureGroups)
#     if test != None:
#         print "Building train examples"
#         testLabels, testFeatures = buildExamples(test, limit, limitTerms, featureGroups)
#     optimize(trainFeatures, develFeatures, trainLabels, develLabels, verbose=3, n_jobs = -1, scoring = "f1_micro", cvJobs=1)

def run(dataPath, outDir=None, actions=None, featureGroups=None, limit=None, numTerms=100, useTestSet=False, clear=False):
    if clear and os.path.exists(outDir):
        print "Removing output directory", outDir
        shutil.rmtree(outDir)
    if not os.path.exists(outDir):
        print "Making output directory", outDir
        os.makedirs(outDir)
    Stream.openLog(os.path.join(options.output, "log.txt"))
    #loadUniprotSimilarity(os.path.join(options.dataPath, "Uniprot", "similar.txt"), proteins)
    terms = loadGOTerms(os.path.join(options.dataPath, "GO", "go_terms.tsv"))
    
    picklePath = os.path.join(outDir, "examples.pickle")
    if actions == None or "build" in actions:
        print "==========", "Building Examples", "=========="
        proteins = defaultdict(lambda: dict())
        loadSequences(os.path.join(options.dataPath, "Swiss_Prot", "Swissprot_sequence.tsv.gz"), proteins)
        termCounts = loadAnnotations(os.path.join(options.dataPath, "Swiss_Prot", "Swissprot_propagated.tsv.gz"), proteins)
        print "Proteins:", len(proteins)
        print "Unique terms:", len(termCounts)
        topTerms = getTopTerms(termCounts, numTerms)
        print "Using", len(topTerms), "most common GO terms"
        #print "Most common terms:", topTerms
        #print proteins["14310_ARATH"]
        loadSplit(os.path.join(options.dataPath, "Swiss_Prot"), proteins)
        #divided = splitProteins(proteins)
        examples = buildExamples(proteins, dataPath, limit, limitTerms=set([x[0] for x in topTerms]), featureGroups=featureGroups)
        saveFeatureNames(examples["feature_names"], os.path.join(outDir, "features.tsv"))
        print "Pickling examples to", picklePath
        with open(picklePath, "wb") as pickeFile:
            pickle.dump(examples, pickeFile)
    else:
        print "Loading examples from", picklePath
        with open(picklePath, "rb") as pickeFile:
            examples = pickle.load(pickeFile)
    if actions == None or "classify" in actions:
        print "==========", "Training Classifier", "=========="
        best = optimize(examples, terms=terms)
        saveResults(best["results"], os.path.join(outDir, "devel-results.tsv"))
    #y, X = buildExamples(proteins, None, set([x[0] for x in topTerms]))
    #print y
    #print X
    #print time.strftime('%X %x %Z')
    #classify(y, X)
    #print time.strftime('%X %x %Z')

if __name__=="__main__":       
    from optparse import OptionParser
    optparser = OptionParser(description="")
    optparser.add_option("-a", "--actions", default=None, help="")
    optparser.add_option("-p", "--dataPath", default=os.path.expanduser("~/data/CAFA3"), help="")
    optparser.add_option("-f", "--features", default="similar", help="")
    optparser.add_option("-l", "--limit", default=None, type=int, help="")
    optparser.add_option("-t", "--terms", default=100, type=int, help="")
    optparser.add_option("-o", "--output", default=None, help="")
    optparser.add_option("--testSet", default=False, action="store_true", help="")
    optparser.add_option("--clear", default=False, action="store_true", help="")
    (options, args) = optparser.parse_args()
    
    if options.actions != None:
        options.actions = options.actions.split(",")
    #proteins = de
    #importProteins(os.path.join(options.dataPath, "Swiss_Prot", "Swissprot_sequence.tsv.gz"))
    run(options.dataPath, actions=options.actions, featureGroups=options.features.split(","), 
        limit=options.limit, numTerms=options.terms, useTestSet=options.testSet, outDir=options.output,
        clear=options.clear)