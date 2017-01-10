import sys, os
import gzip
from collections import defaultdict
import csv
from sklearn.feature_extraction import DictVectorizer
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.grid_search import GridSearchCV, ParameterGrid
from sklearn.ensemble.forest import ExtraTreesClassifier, RandomForestClassifier
from sklearn.metrics import classification_report, f1_score, precision_score, recall_score
from featureBuilders import *
from utils import Stream
import operator
import time
from sklearn.cross_validation import train_test_split
from sklearn.metrics.ranking import roc_auc_score
import shutil
import cPickle as pickle
from sklearn.multiclass import OneVsRestClassifier

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

def addProtein(proteins, protId, cafaId, sequence, filename):
    assert len(sequence) > 0
    if protId not in proteins:
        proteins[protId]["seq"] = sequence
        proteins[protId]["id"] = protId
        proteins[protId]["cafa_id"] = cafaId
        proteins[protId]["file"] = [filename]
    else:
        if proteins[protId]["seq"] != sequence:
            print "WARNING, sequence mismatch for", protId, cafaId, ":", (proteins[protId], (protId, cafaId, sequence, filename))
        assert proteins[protId]["id"] == protId, (proteins[protId], (protId, cafaId, sequence, filename))
        proteins[protId]["file"] += [filename]
        if cafaId != None:
            if proteins[protId]["cafa_id"] != None:
                print "WARNING, duplicate CAFA target", protId, cafaId, ":", (proteins[protId], (protId, cafaId, sequence, filename))
            proteins[protId]["cafa_id"] = cafaId
         
def loadFASTA(inPath, proteins, cafaHeader=False):
    print "Loading sequences from", inPath
    filename = os.path.basename(inPath)
    with openAny(inPath, "rt") as f:
        protId = None
        cafaId = None
        sequence = ""
        for line in f:
            if line.startswith(">"):
                # Add already read protein
                if protId != None:
                    addProtein(proteins, protId, cafaId, sequence, filename)
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
            addProtein(proteins, protId, cafaId, sequence, filename)               
            #print seq.id, seq.seq

def loadSplit(inPath, proteins):
    for dataset in ("train", "devel", "test"):
        filePath = os.path.join(inPath, dataset + ".txt.gz")
        assert os.path.exists(filePath), filePath
        with gzip.open(filePath, "rt") as f:
            for line in f:
                protId = line.strip()
                assert protId in proteins
                proteins[protId]["origSet"] = dataset

def defineSets(proteins, cafaTargets):
    counts = defaultdict(int)
    for protein in proteins.values():
        cafaSet = ["cafa"] if protein["cafa_id"] != None else []
        origSet = [protein["origSet"]] if protein.get("origSet") != None else []
        if len(cafaSet) > 0:
            if cafaTargets == "overlap":
                protein["sets"] = cafaSet + origSet
            elif cafaTargets == "separate":
                protein["sets"] = cafaSet
            else:
                raise Exception("CAFA targets were loaded with mode '" + cafaTargets + "'")
        else:
            protein["sets"] = origSet
        category = ",".join(origSet + cafaSet) + "=>" + ",".join(protein["sets"])
        counts[category] += 1
    print "Defined sets:", dict(counts)

def saveFeatureNames(names, outPath):
    print "Saving feature names to", outPath
    with open(outPath, "wt") as f:
        f.write("index\tname\n")
        for i in range(len(names)):
            f.write(str(i) + "\t" + names[i] + "\n")
    
def vectorizeExamples(examples):
    mlb = MultiLabelBinarizer()
    dv = DictVectorizer(sparse=True)
    examples["features"] = dv.fit_transform(examples["features"])
    examples["feature_names"] = dv.feature_names_
    examples["labels"] = mlb.fit_transform(examples["labels"])
    examples["label_names"] = mlb.classes_
    print "Vectorized", len(examples["labels"]), "examples with", len(examples["feature_names"]), "unique features and", len(examples["label_names"]), "unique labels"
        
# def splitProteins(proteins):
#     datasets = {"devel":[], "train":[], "test":[]}
#     for protId in sorted(proteins.keys()):
#         datasets[proteins[protId]["set"]].append(proteins[protId])
#     print "Divided sets", [(x, len(datasets[x])) for x in sorted(datasets.keys())]
#     return datasets

def buildExamples(proteins, dataPath, limit=None, limitTerms=None, featureGroups=None):
    print "Building examples"
    examples = {"labels":[], "features":[], "ids":[], "sets":[], "label_names":[], "label_size":{}}
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
    print "Building features, feature groups = ", featureGroups
    if featureGroups == None or "taxonomy" in featureGroups:
        builder = TaxonomyFeatureBuilder(os.path.join(dataPath, "Taxonomy"))
        builder.build(protObjs)
    if featureGroups == None or "similar" in featureGroups:
        builder = UniprotFeatureBuilder(os.path.join(dataPath, "Uniprot", "similar.txt"))
        builder.build(protObjs)
    if featureGroups == None or "blast" in featureGroups:
        builder = BlastFeatureBuilder(os.path.join(dataPath, "blastp_result_features"))
        builder.build(protObjs)
    builder = None
    examples["features"] = [x["features"] for x in protObjs]
    for protObj in protObjs:
        del protObj["features"]
    # Prepare the examples
    print "Built", len(examples["labels"]), "examples" # with", len(examples["feature_names"]), "unique features"
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

def metricsToString(result, style="%.3f"):
    return "a/f|p/r|tp/fp/tn/fn = " + style % result["auc"] + "/" + style % result["fscore"] + "|" + style % result["precision"] + "/" + style % result["recall"] \
        + "|" + "/".join([str(result.get(x, "-")) for x in ("tp", "fp", "tn", "fn")])

def getResultsString(results, maxNumber=None, skipIds=None):
    count = 0
    s = ""
    for result in sorted(results.values(), key=lambda x: x["auc"], reverse=True):
        if skipIds != None and result["id"] in skipIds:
            continue
        s += metricsToString(result) + " " + str([result.get("id"), result.get("ns"), result.get("label_size"), result.get("name")]) + "\n"
        count += 1
        if count > maxNumber:
            break
    return s

def saveResults(results, outPath):
    print "Writing results to", outPath
    with open(outPath, "wt") as f:
        dw = csv.DictWriter(f, ["auc", "fscore", "precision", "recall", "tp", "fp", "tn", "fn", "id", "label_size", "ns", "name"], delimiter='\t')
        dw.writeheader()
        dw.writerow(results["average"])
        results = [x for x in results.values() if x["id"] != "average"]
        dw.writerows(sorted(results, key=lambda x: x["auc"], reverse=True))

def getMatch(gold, predicted):
    if gold == predicted:
        return "tp" if (gold == 1) else "tn"
    else:
        return "fn" if (gold == 1) else "fp"

def savePredictions(exampleIds, labels, predicted, label_names, outPath):
    print "Writing predictions to", outPath
    lengths = [len(x) for x in (exampleIds, labels, predicted)]
    assert len(set(lengths)) == 1, lengths
    label_indices = range(len(label_names))
    rows = []
    for exampleId, gold, pred in zip(exampleIds, labels, predicted):
        for i in label_indices:
            if gold[i] == 1 or pred[i] == 1:
                row = {"id":exampleId, "label":label_names[i], "gold":gold[i], "predicted":int(pred[i])}
                row["match"] = getMatch(gold[i], pred[i])
                rows.append(row)
    with open(outPath, "wt") as f:
        dw = csv.DictWriter(f, ["id", "label", "gold", "predicted", "match"], delimiter='\t')
        dw.writeheader()
        dw.writerows(rows)

def importNamed(name):
    asName = name.rsplit(".", 1)[-1]
    imported = False
    attempts = ["from sklearn." + name.rsplit(".", 1)[0] + " import " + asName,
                "from " + name.rsplit(".", 1)[0] + " import " + asName,
                "import " + name + " as " + asName]
    for attempt in attempts:
        try:
            print "Importing '" + attempt + "', ",
            exec attempt
            imported = True
            print "OK"
            break;
        except ImportError:
            print "failed"
    if not imported:
        raise Exception("Could not import '" + name + "'")
    return eval(asName)

# def parseOptions(string):
#     string = string.strip()
#     if not string.startswith("{"):
#         string = "{" + string
#     if not string.endswith("}"):
#         string = "}" + string
#     print "Parsing options string:", string
#     return eval(string)

def evaluate(labels, predicted, label_names, label_size=None, terms=None):
    print "Evaluating the predictions"
    results = {}
    # Get the average result
    results["average"] = {"id":"average", "ns":None, "name":None, "auc":0, "tp":None, "fp":None, "fn":None, "tn":None}
    try:
        results["average"]["auc"] = roc_auc_score(labels, predicted, average="micro")
    except ValueError as e:
        print e
    results["average"]["fscore"] = f1_score(labels, predicted, average="micro")
    results["average"]["precision"] = precision_score(labels, predicted, average="micro")
    results["average"]["recall"] = recall_score(labels, predicted, average="micro")
    # Get results per label
    try:
        aucs = roc_auc_score(labels, predicted, average=None)
    except ValueError as e:
        print e
        aucs = [0] * len(label_names)
    fscores = f1_score(labels, predicted, average=None)
    precisions = precision_score(labels, predicted, average=None)
    recalls = recall_score(labels, predicted, average=None)
    lengths = [len(x) for x in (aucs, fscores, precisions, recalls, label_names)]
    assert len(set(lengths)) == 1, lengths
    for auc, fscore, precision, recall, label_name in zip(aucs, fscores, precisions, recalls, label_names):
        assert label_name not in results
        result = {"id":label_name, "ns":None, "name":None, "auc":auc, "precision":precision, "recall":recall, "fscore":fscore, "tp":0, "fp":0, "fn":0, "tn":0}
        results[label_name] = result
        if label_size != None:
            result["label_size"] = label_size[label_name]
        if terms != None and label_name in terms:
            term = terms[label_name]
            result["ns"] = term["ns"]
            result["name"] = term["name"]
    # Calculate instances
    stats = {x:{"tp":0, "fp":0, "fn":0, "tn":0} for x in label_names}
    label_indices = range(len(label_names))
    for gold, pred in zip(labels, predicted):
        for i in label_indices:
            stats[label_names[i]][getMatch(gold[i], pred[i])] += 1
    for key in stats:
        results[key].update(stats[key])
    return results

def optimize(classifier, classifierArgs, examples, cvJobs=1, terms=None, useOneVsRest=False, outDir=None):
    #grid = ParameterGrid({"n_estimators":[10], "n_jobs":[n_jobs], "verbose":[verbose]}) #{"n_estimators":[1,2,10,50,100]}
    #XTrainAndDevel, XTest, yTrainAndDevel, yTest = train_test_split(X, y, test_size=0.2, random_state=0)
    #XTrain, XDevel, yTrain, yDevel = train_test_split(XTrainAndDevel, yTrainAndDevel, test_size=0.2, random_state=0)
    sets = examples["sets"]
    trainIndices = [i for i in range(len(sets)) if "train" in sets[i]]
    develIndices = [i for i in range(len(sets)) if "devel" in sets[i]]
    trainFeatures = examples["features"][trainIndices]
    develFeatures = examples["features"][develIndices]
    trainLabels = examples["labels"][trainIndices]
    develLabels = examples["labels"][develIndices]
    develIds = [examples["ids"][i] for i in range(len(sets)) if sets[i] == "devel"]
    print "Optimizing, train / devel = ", trainFeatures.shape[0], "/", develFeatures.shape[0]
    if useOneVsRest:
        print "Using OneVsRestClassifier"
    best = None
    print "Parameter grid search", time.strftime('%X %x %Z')
    Cls = importNamed(classifier)
    #grid = parseOptions(classifierArgs)
    for args in ParameterGrid(classifierArgs):
        print "Learning with args", args
        cls = Cls(**args)
        if useOneVsRest:
            cls = OneVsRestClassifier(cls)
        cls.fit(trainFeatures, trainLabels)
        print "Predicting the devel set"
        predicted = cls.predict(develFeatures)
        #score = roc_auc_score(develLabels, predicted, average="micro")
        #scores = roc_auc_score(develLabels, predicted, average=None)
        #print "Average =", score
        #results = getResults(examples, scores, terms)
        results = evaluate(develLabels, predicted, examples["label_names"], examples["label_size"], terms)
        print "Average:", metricsToString(results["average"])
        print getResultsString(results, 20, ["average"])
        if best == None or results["average"]["auc"] > best["results"]["average"]["auc"]:
            best = {"results":results, "args":args, "predicted":predicted, "gold":develLabels}
            if hasattr(cls, "feature_importances_"):
                best["feature_importances"] = cls.feature_importances_
        print time.strftime('%X %x %Z')
    if outDir != None:
        saveResults(best["results"], os.path.join(outDir, "devel-results.tsv"))
        savePredictions(develIds, develLabels, best["predicted"], examples["label_names"], os.path.join(outDir, "devel-predictions.tsv"))
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

def run(dataPath, outDir=None, actions=None, featureGroups=None, classifier=None, classifierArgs=None, useOneVsRest=False, limit=None, numTerms=100, useTestSet=False, clear=False, cafaTargets="skip"):
    if clear and os.path.exists(outDir):
        print "Removing output directory", outDir
        shutil.rmtree(outDir)
    if not os.path.exists(outDir):
        print "Making output directory", outDir
        os.makedirs(outDir)
    Stream.openLog(os.path.join(options.output, "log.txt"))
    if actions != None:
        for action in actions:
            assert action in ("build", "classify")
    assert cafaTargets in ("skip", "overlap", "separate")
    #loadUniprotSimilarity(os.path.join(options.dataPath, "Uniprot", "similar.txt"), proteins)
    terms = loadGOTerms(os.path.join(options.dataPath, "GO", "go_terms.tsv"))
    
    picklePath = os.path.join(outDir, "examples.pickle.gz")
    examples = None
    if actions == None or "build" in actions:
        print "==========", "Building Examples", "=========="
        proteins = defaultdict(lambda: dict())
        loadFASTA(os.path.join(options.dataPath, "Swiss_Prot", "Swissprot_sequence.tsv.gz"), proteins)
        if cafaTargets != "skip":
            print "Loading CAFA3 targets"
            loadFASTA(os.path.join(options.dataPath, "CAFA3_targets", "Target_files", "target.all.fasta"), proteins, True)
        print "Proteins:", len(proteins)
        termCounts = loadAnnotations(os.path.join(options.dataPath, "Swiss_Prot", "Swissprot_propagated.tsv.gz"), proteins)
        print "Unique terms:", len(termCounts)
        topTerms = getTopTerms(termCounts, numTerms)
        print "Using", len(topTerms), "most common GO terms"
        #print "Most common terms:", topTerms
        #print proteins["14310_ARATH"]
        loadSplit(os.path.join(options.dataPath, "Swiss_Prot"), proteins)
        defineSets(proteins, cafaTargets)
        #divided = splitProteins(proteins)
        examples = buildExamples(proteins, dataPath, limit, limitTerms=set([x[0] for x in topTerms]), featureGroups=featureGroups)
        print "Pickling examples to", picklePath
        with gzip.open(picklePath, "wb") as pickleFile:
            pickle.dump(examples, pickleFile)
    if actions == None or "classify" in actions:
        print "==========", "Training Classifier", "=========="
        if examples == None:
            print "Loading examples from", picklePath
            with gzip.open(picklePath, "rb") as pickleFile:
                examples = pickle.load(pickleFile)
        vectorizeExamples(examples)
        if not os.path.exists(os.path.join(outDir, "features.tsv")):
            saveFeatureNames(examples["feature_names"], os.path.join(outDir, "features.tsv"))
        best = optimize(classifier, classifierArgs, examples, terms=terms, useOneVsRest=useOneVsRest, outDir=outDir)
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
    optparser.add_option("-f", "--features", default="blast", help="")
    optparser.add_option("-l", "--limit", default=None, type=int, help="")
    optparser.add_option("-t", "--terms", default=100, type=int, help="")
    optparser.add_option("-o", "--output", default=None, help="")
    optparser.add_option('-c','--classifier', help='', default="ensemble.RandomForestClassifier")
    optparser.add_option('-r','--args', help='', default="{'n_estimators':[10], 'n_jobs':[1], 'verbose':[3]}")
    optparser.add_option("--onevsrest", default=False, action="store_true", help="")
    optparser.add_option("--testSet", default=False, action="store_true", help="")
    optparser.add_option("--clear", default=False, action="store_true", help="")
    optparser.add_option("--targets", default="skip", help="skip, overlap or separate")
    (options, args) = optparser.parse_args()
    
    if options.actions != None:
        options.actions = options.actions.split(",")
    options.args = eval(options.args)
    #proteins = de
    #importProteins(os.path.join(options.dataPath, "Swiss_Prot", "Swissprot_sequence.tsv.gz"))
    run(options.dataPath, actions=options.actions, featureGroups=options.features.split(","), 
        limit=options.limit, numTerms=options.terms, useTestSet=options.testSet, outDir=options.output,
        clear=options.clear, classifier=options.classifier, classifierArgs=options.args, 
        useOneVsRest=options.onevsrest, cafaTargets=options.targets)