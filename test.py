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
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2
from sklearn.feature_selection.variance_threshold import VarianceThreshold
from sklearn.multioutput import MultiOutputClassifier, _fit_estimator
from sklearn.utils.validation import check_X_y, has_fit_parameter
from sklearn.externals.joblib.parallel import Parallel, delayed
import numpy as np
try:
    import ujson as json
except ImportError:
    import json
import statistics

class MyMultiOutputClassifier(MultiOutputClassifier):
    def fit(self, X, y, sample_weight=None):
        if not hasattr(self.estimator, "fit"):
            raise ValueError("The base estimator should implement a fit method")
    
        X, y = check_X_y(X, y,
                         multi_output=True,
                         accept_sparse=True)
    
        if y.ndim == 1:
            raise ValueError("y must have at least two dimensions for "
                             "multi target regression but has only one.")
    
        if (sample_weight is not None and
                not has_fit_parameter(self.estimator, 'sample_weight')):
            raise ValueError("Underlying regressor does not support"
                             " sample weights.")
    
        self.estimators_ = Parallel(n_jobs=self.n_jobs, verbose=3)(delayed(_fit_estimator)(
            self.estimator, X, y[:, i], sample_weight) for i in range(y.shape[1]))
        return self

#MultiOutputClassifier.fit = newFit

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
        
# def splitProteins(proteins):
#     datasets = {"devel":[], "train":[], "test":[]}
#     for protId in sorted(proteins.keys()):
#         datasets[proteins[protId]["set"]].append(proteins[protId])
#     print "Divided sets", [(x, len(datasets[x])) for x in sorted(datasets.keys())]
#     return datasets

def getFeatureGroups(groups=None):
    if groups == None:
        groups = ["all"]
    groups = list(set(groups))
    if "all" in groups:
        groups.remove("all")
        groups += set(["taxonomy", "blast", "delta", "interpro", "predgpi", "nucpred"])
    removed = [x for x in groups if x.startswith("-")]
    groups = [x for x in groups if not x.startswith("-")]
    for group in removed:
        groups.remove(group.strip("-"))
    return groups

def buildExamples(proteins, dataPath, limit=None, limitTerms=None, featureGroups=None):
    print "Building examples"
    examples = {"labels":[], "features":[], "ids":[], "cafa_ids":[], "sets":[], "label_names":[], "label_size":{}}
    protIds = sorted(proteins.keys())
    if limit:
        protIds = protIds[0:limit]
    protObjs = [proteins[key] for key in protIds]
    for protein in protObjs:
        # Initialize features
        protein["features"] = {"DUMMY:dummy":1}
        # Build labels
        labels = protein["terms"].keys()
        if limitTerms:
            labels = [x for x in labels if x in limitTerms]
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
    # Build features
    featureGroups = getFeatureGroups(featureGroups)
    for group in featureGroups:
        if group not in ("taxonomy", "similar", "blast", "delta", "interpro", "predgpi", "nucpred"):
            raise Exception("Unknown feature group '" + str(group) + "'")
    print "Building features, feature groups =", featureGroups
    if featureGroups == None or "taxonomy" in featureGroups:
        builder = TaxonomyFeatureBuilder([os.path.join(dataPath, "Taxonomy")])
        builder.build(protObjs)
    if featureGroups == None or "similar" in featureGroups:
        builder = UniprotFeatureBuilder(os.path.join(dataPath, "Uniprot", "similar.txt"))
        builder.build(protObjs)
    if featureGroups == None or "blast" in featureGroups:
        builder = BlastFeatureBuilder([os.path.join(dataPath, "temp_blastp_result_features"), os.path.join(dataPath, "blastp_result_features")])
        builder.build(protObjs)
    if featureGroups == None or "delta" in featureGroups:
        builder = BlastFeatureBuilder([os.path.join(dataPath, "temp_deltablast_result_features"), os.path.join(dataPath, "deltablast_result_features")], tag="DELTA")
        builder.build(protObjs)
    if featureGroups == None or "interpro" in featureGroups:
        builder = InterProScanFeatureBuilder([os.path.join(dataPath, "temp_interproscan_result_features"), os.path.join(dataPath, "interproscan_result_features")])
        builder.build(protObjs)
    if featureGroups == None or "predgpi" in featureGroups:
        builder = GPIAnchoringFeatureBuilder([os.path.join(dataPath, "predGPI")])
        builder.build(protObjs)
    if featureGroups == None or "nucpred" in featureGroups:
        builder = NucPredFeatureBuilder([os.path.join(dataPath, "nucPred")])
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

def learn(Cls, args, examples, trainSets, testSets, useMultiOutputClassifier, terms, cls=None, averageOnly=False):
    args = args.copy()
    print "Learning with args", args
    if cls == None:
        print "Initializing classifier", Cls.__name__, "with arguments", args
        cls = Cls(**args)
        if useMultiOutputClassifier:
            print "Using MultiOutputClassifier"
            cls = MyMultiOutputClassifier(cls) # OneVsRestClassifier(cls)
    else:
        print "Using existing classifier"
    print "Train / test = ", trainSets, "/", testSets
    sets = examples["sets"]
    trainIndices = [i for i in range(len(sets)) if any(x in trainSets for x in sets[i])]
    testIndices = [i for i in range(len(sets)) if any(x in testSets for x in sets[i])]
    trainFeatures = examples["features"][trainIndices]
    testFeatures = examples["features"][testIndices]
    trainLabels = examples["labels"][trainIndices]
    testLabels = examples["labels"][testIndices]
    testIds = [examples["ids"][i] for i in range(len(sets)) if any(x in testSets for x in sets[i])]
    testCafaIds = [examples["cafa_ids"][i] for i in range(len(sets)) if any(x in testSets for x in sets[i])]
    print "Training, train / test = ", trainFeatures.shape[0], "/", testFeatures.shape[0]
    cls.fit(trainFeatures, trainLabels)
    print "Predicting"
    predicted = cls.predict(testFeatures)
    probabilities = None
    if hasattr(cls, "predict_proba"):
        print "Predicting probabilities"
        probabilities = cls.predict_proba(testFeatures)
    results = evaluate(testLabels, predicted, examples["label_names"], examples["label_size"], terms, averageOnly=averageOnly)
    print "Average:", metricsToString(results["average"])
    if not averageOnly:
        print getResultsString(results, 20, ["average"])
    #if predictionsPath != None:
    #    predictionsPath(testIds, testLabels, predicted, examples["label_names"], predictionsPath)
    data = {"results":results, "args":args, "predicted":predicted, "gold":testLabels, "ids":testIds, "cafa_ids":testCafaIds, "probabilities":probabilities}
    if hasattr(cls, "feature_importances_"):
        data["feature_importances"] = cls.feature_importances_        
    return cls, data

def resultIsBetter(original, new, key="average"):
    if new[key]["fscore"] != original[key]["fscore"]:
        return new[key]["fscore"] > original[key]["fscore"]
    else:
        return new[key]["auc"] > original[key]["auc"]

def warmStartGrid(Cls, classifierArgs, examples, terms):
    print "Using warm start parameter grid search"
    for key in classifierArgs:
        if key != "n_estimators" and len(classifierArgs[key]) > 1:
            raise Exception("Multiple classifier argument values defined for argument '" + str(key) + "'")
    args = classifierArgs.copy()
    numEstimatorList = sorted(args["n_estimators"])
    args["n_estimators"] = [numEstimatorList[0]]
    args = {x:args[x][0] for x in args.keys()}
    cls = None
    best = None
    performances = []
    for n in numEstimatorList:
        if cls != None:
            args["n_estimators"] = n
            cls.n_estimators = n
            print "cls.n_estimators = ", cls.n_estimators
        cls, data = learn(Cls, args, examples, ["train"], ["devel"], False, terms, cls=cls, averageOnly=True)
        performances.append({x:data["results"]["average"][x] for x in ("auc", "fscore", "precision", "recall")})
        performances[-1]["n"] = n
        if best == None or resultIsBetter(best["results"], data["results"]):
            best = data
        else: # Release the not-best results
            data = None
    print "Warm start parameter grid search complete"
    for performance in performances:
        print performance["n"], "\t", metricsToString(performance)
    print "Full evaluation for the best results"
    best["results"] = evaluate(best["gold"], best["predicted"], examples["label_names"], examples["label_size"], terms)
    return best

def optimize(classifier, classifierArgs, examples, cvJobs=1, terms=None, useMultiOutputClassifier=False, outDir=None, negatives=False, useTestSet=False, useCAFASet=False):
    best = None
    print "Parameter grid search"
    Cls = importNamed(classifier)
    if classifierArgs.get("warm_start") == [True]:
        best = warmStartGrid(Cls, classifierArgs, examples, terms)
    else:
        for args in ParameterGrid(classifierArgs):
            _, data = learn(Cls, args, examples, ["train"], ["devel"], useMultiOutputClassifier, terms)
            if best == None or resultIsBetter(best["results"], data["results"]):
                best = data #{"results":results, "args":args, "predicted":predicted, "gold":develLabels}
            else: # Release the not-best results
                data = None
    print "Best classifier arguments:", best["args"]
    print "Best development set results:", metricsToString(best["results"]["average"])
    print getResultsString(best["results"], 20, ["average"])
    if outDir != None:
        saveResults(best, os.path.join(outDir, "devel"), examples["label_names"], negatives=negatives)
    if useTestSet:
        print "Classifying the test set"
        _, data = learn(Cls, best["args"], examples, ["train", "devel"], ["test"], useMultiOutputClassifier, terms)
        if outDir != None:
            saveResults(data, os.path.join(outDir, "test"), examples["label_names"], negatives=negatives)
    if useCAFASet:
        print "Classifying the CAFA targets"
        _, data = learn(Cls, best["args"], examples, ["train", "devel", "test"], ["cafa"], useMultiOutputClassifier, terms)
        if outDir != None:
            saveResults(data, os.path.join(outDir, "cafa"), examples["label_names"], negatives=negatives)

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

def run(dataPath, outDir=None, actions=None, featureGroups=None, classifier=None, classifierArgs=None, useMultiOutputClassifier=False, limit=None, numTerms=100, useTestSet=False, clear=False, cafaTargets="skip", negatives=False):
    if clear and os.path.exists(outDir):
        print "Removing output directory", outDir
        shutil.rmtree(outDir)
    if not os.path.exists(outDir):
        print "Making output directory", outDir
        os.makedirs(outDir)
    Stream.openLog(os.path.join(options.output, "log.txt"))
    if actions != None:
        for action in actions:
            assert action in ("build", "classify", "statistics")
    assert cafaTargets in ("skip", "overlap", "separate", "external")
    #loadUniprotSimilarity(os.path.join(options.dataPath, "Uniprot", "similar.txt"), proteins)
    terms = loadGOTerms(os.path.join(options.dataPath, "GO", "go_terms.tsv"))
    
    exampleFilePath = os.path.join(outDir, "examples.json.gz")
    examples = None
    if actions == None or "build" in actions:
        print "==========", "Building Examples", "=========="
        proteins = defaultdict(lambda: dict())
        loadFASTA(os.path.join(options.dataPath, "Swiss_Prot", "Swissprot_sequence.tsv.gz"), proteins)
        if cafaTargets != "skip":
            print "Loading CAFA3 targets"
            loadFASTA(os.path.join(options.dataPath, "CAFA3_targets", "Target_files", "target.all.fasta"), proteins, True)
        print "Proteins:", len(proteins)
        termCounts = loadAnnotations(os.path.join(options.dataPath, "data", "Swissprot_propagated.tsv.gz"), proteins)
        print "Unique terms:", len(termCounts)
        topTerms = getTopTerms(termCounts, numTerms)
        print "Using", len(topTerms), "most common GO terms"
        #print "Most common terms:", topTerms
        #print proteins["14310_ARATH"]
        loadSplit(os.path.join(options.dataPath, "data"), proteins)
        defineSets(proteins, cafaTargets)
        #divided = splitProteins(proteins)
        examples = buildExamples(proteins, dataPath, limit, limitTerms=set([x[0] for x in topTerms]), featureGroups=featureGroups)
        print "Saving examples to", exampleFilePath
        with gzip.open(exampleFilePath, "wt") as pickleFile:
            json.dump(examples, pickleFile, indent=2) #pickle.dump(examples, pickleFile)
        vectorizeExamples(examples, featureGroups)
    elif len([x for x in actions if x != "build"]) > 0:
        print "==========", "Loading Examples", "=========="
        if examples == None:
            print "Loading examples from", exampleFilePath
            with gzip.open(exampleFilePath, "rt") as pickleFile:
                examples = json.load(pickleFile) #pickle.load(pickleFile)
        vectorizeExamples(examples, featureGroups)
    if actions == None or "classify" in actions:
        print "==========", "Training Classifier", "=========="
        if not os.path.exists(os.path.join(outDir, "features.tsv")):
            saveFeatureNames(examples["feature_names"], os.path.join(outDir, "features.tsv"))
        optimize(classifier, classifierArgs, examples, terms=terms, 
                 useMultiOutputClassifier=useMultiOutputClassifier, outDir=outDir, negatives=negatives,
                 useTestSet=useTestSet, useCAFASet=(cafaTargets != "skip"))
    if actions == None or "statistics" in actions:
        print "==========", "Calculating Statistics", "=========="
        statistics.makeStatistics(examples, outDir)
    #y, X = buildExamples(proteins, None, set([x[0] for x in topTerms]))
    #print y
    #print X
    #print time.strftime('%X %x %Z')
    #classify(y, X)
    #print time.strftime('%X %x %Z')

if __name__=="__main__":       
    from optparse import OptionParser
    optparser = OptionParser(description="")
    optparser.add_option("-a", "--actions", default=None, help="One or more of 'build', 'classify' or 'statistics'")
    optparser.add_option("-p", "--dataPath", default=os.path.expanduser("~/data/CAFA3/data"), help="The main directory for the data files")
    optparser.add_option("-f", "--features", default="all", help="Comma-separated list of feature group names. Use 'all' for all feature groups and '-name' to remove groups.")
    optparser.add_option("-l", "--limit", default=None, type=int, help="Limit the number of proteins to read.")
    optparser.add_option("-t", "--terms", default=100, type=int, help="The number of top most common GO terms to use as labels")
    optparser.add_option("-o", "--output", default=None, help="The output directory")
    optparser.add_option('-c','--classifier', help='', default="ensemble.RandomForestClassifier")
    optparser.add_option('-r','--args', help='', default="{'random_state':[1], 'n_estimators':[10], 'n_jobs':[1], 'verbose':[3]}")
    optparser.add_option("--multioutputclassifier", default=False, action="store_true", help="Use the MultiOutputClassifier to train a separate classifier for each label")
    optparser.add_option("--testSet", default=False, action="store_true", help="Classify the test set")
    optparser.add_option("--clear", default=False, action="store_true", help="Remove the output directory if it already exists")
    optparser.add_option("--targets", default="skip", help="How to include the CAFA target proteins, one of 'skip', 'overlap' or 'separate'")
    optparser.add_option("--negatives", default=False, action="store_true", help="Write negative predictions in the result files")
    (options, args) = optparser.parse_args()
    
    if options.actions != None:
        options.actions = options.actions.split(",")
    options.args = eval(options.args)
    #proteins = de
    #importProteins(os.path.join(options.dataPath, "Swiss_Prot", "Swissprot_sequence.tsv.gz"))
    run(options.dataPath, actions=options.actions, featureGroups=options.features.split(","), 
        limit=options.limit, numTerms=options.terms, useTestSet=options.testSet, outDir=options.output,
        clear=options.clear, classifier=options.classifier, classifierArgs=options.args, 
        useMultiOutputClassifier=options.multioutputclassifier, cafaTargets=options.targets,
        negatives=options.negatives)