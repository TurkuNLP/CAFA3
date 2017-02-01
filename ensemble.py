from sklearn.feature_extraction import DictVectorizer
from collections import defaultdict
import evaluateFile
import loading
import evaluation
import sys, os
import classification
import shutil
from utils import Stream
from sklearn.grid_search import GridSearchCV
import itertools

def clearKey(proteins, key):
    for protId in proteins:
        protein = proteins[protId]
        if key in protein:
            del protein[key]

def getCombinations(items):
    combinations = []
    for i in xrange(1, len(items) + 1):
        els = [list(x) for x in itertools.combinations(items, i)]
        combinations.extend(els)
    return combinations

def combinePred(proteins, predKeys, combKey, mode="AND", limitToSets=None):
    assert mode in ("AND", "OR", "SINGLE")
    if mode == "SINGLE":
        assert len(predKeys) == 1
    counts = defaultdict(int)
    for protId in proteins:
        protein = proteins[protId]
        if limitToSets != None and not any(x in limitToSets for x in protein["sets"]):
            counts["out-of-sets"] += 1
            continue
        counts["proteins"] += 1
        protein[combKey] = {}
        if mode == "OR":
            for key in predKeys:
                if key in protein:
                    counts["predictions-mode-" + mode] += 1
                    protein[combKey].update(protein[key])
                else:
                    counts["no-prediction-for-" + key] += 1
        elif mode == "AND":
            missing = False
            predLabelSets = []
            for key in predKeys:
                if key not in protein:
                    missing = True
                    counts["no-prediction-for-" + key] += 1
                else:
                    predLabelSets.append(set(protein[key].keys()))
            if not missing:
                counts["predictions-mode-" + mode] += 1
                #pred1 = set(protein[key1].keys())
                #pred2 = set(protein[key2].keys())
                protein[combKey] = {x:1 for x in set.intersection(predLabelSets)} #{x:1 for x in pred1.intersection(pred2)}
            else:
                protein[combKey] = {}
        else:
            key = predKeys[0] #key1 if mode == "ONLY1" else key2
            if key not in protein:
                counts["no-prediction-for-" + key] += 1
            else:
                counts["predictions-mode-" + mode] += 1
                protein[combKey] = protein[key]
                
    print "Combined predictions, mode =", mode, "counts =", dict(counts)

def buildExamples(proteins, key1, key2, limitToSets=None, limitTerms=None, outDir=None):
    counts = defaultdict(int)
    empty = {}
    examples = {"classes":[], "features":[], "sets":[], "proteins":[], "labels":[]}
    for protId in proteins:
        protein = proteins[protId]
        if limitToSets != None and not any(x in limitToSets for x in protein["sets"]):
            counts["out-of-sets"] += 1
            continue
        if key1 not in protein:
            counts["no-prediction-for-" + key1] += 1
        if key2 not in protein:
            counts["no-prediction-for-" + key2] += 1
        pred1 = protein.get(key1, empty)
        pred2 = protein.get(key2, empty)
        conf1 = protein.get(key1 + "_conf", empty)
        conf2 = protein.get(key2 + "_conf", empty)
        protSets = protein.get("sets")
        goldLabels = protein["terms"].keys()
        if limitTerms:
            goldLabels = [x for x in goldLabels if x in limitTerms]
        allLabels = sorted(set(pred1.keys() + pred2.keys())) # + goldLabels))
        goldLabels = set(goldLabels)
        #for label in goldLabels:
        #    if label not in examples["label_size"]:
        #        examples["label_size"][label] = 0
        #    examples["label_size"][label] += 1
        for label in allLabels:
            features = {} #{label:1}
            features[label] = 1
            for key, pred, conf in ((key1, pred1, conf1), (key2, pred2, conf2)):
                if label in pred:
                    #assert label in conf, (key, pred, conf, counts)
                    features["pos:" + key] = 1
                    features["conf:" + key] = conf[label]
                else:
                    features["neg:" + key] = 1
            cls = 1 if label in goldLabels else 0
            counts["examples"] += 1
            counts["pos" if cls == 1 else "neg"] += 1
            examples["classes"].append(cls)
            examples["features"].append(features)
            examples["sets"].append(protSets)
            examples["proteins"].append(protein)
            examples["labels"].append(label)
    print "Built examples,", dict(counts)
    dv = DictVectorizer(sparse=True)
    examples["features"] = dv.fit_transform(examples["features"])
    examples["feature_names"] = dv.feature_names_
    print "Vectorized the examples, unique features =", len(examples["feature_names"])
    if outDir != None:
        loading.saveFeatureNames(examples["feature_names"], os.path.join(outDir, "features.tsv"))
    return examples

def getSubset(examples, setNames):
    subset = {}
    counts = {}
    sets = examples["sets"]
    indices = [i for i in range(len(sets)) if any(x in setNames for x in sets[i])]
    subset["features"] = examples["features"][indices]
    counts["features"] = subset["features"].shape[0] 
    for key in ("classes", "sets", "proteins", "labels"):
        subset[key] = [examples[key][i] for i in indices]
        counts[key] = len(subset[key])
    print "Generated example subset for sets", setNames, "with", counts
    return subset

def learn(examples, Classifier, classifierArgs, develFolds=10, verbose=3, n_jobs=1, predKey="ensemble", limitTerms=None):
    print "Parameter grid search"
    develExamples = getSubset(examples, ["devel"])
    clf = GridSearchCV(Classifier(), classifierArgs, cv=develFolds, verbose=verbose, n_jobs=n_jobs, scoring="f1_micro")
    clf.fit(develExamples["features"], develExamples["classes"])
    print "Best params", (clf.best_params_, clf.best_score_)
    print "Predicting the test set"
    testExamples = getSubset(examples, ["test"])
    testPredictions = clf.predict(testExamples["features"])
    testProbabilities = clf.predict_proba(testExamples["features"])
    predKeyConf = predKey + "_conf"
    print "Converting binary predictions to labels"
    for prediction, probability, protein, label in zip(testPredictions, testProbabilities, testExamples["proteins"], testExamples["labels"]):
        if predKey not in protein:
            protein[predKey] = {}
            protein[predKeyConf] = {}
        if prediction == 1:
            protein[predKey][label] = 1
            protein[predKeyConf][label] = probability
    print "Evaluating test set ensemble predictions"
    proteins = {x["id"]:x for x in testExamples["proteins"]}
    examples = evaluateFile.makeExamples(proteins, limitTerms=limitTerms, limitToSets=["test"], predKey=predKey)
    loading.vectorizeExamples(examples, None)
    results = evaluation.evaluate(examples["labels"], examples["predictions"], examples, terms=None, averageOnly=True)
    print "Average:", evaluation.metricsToString(results["average"])
    
def combine(dataPath, nnInput, clsInput, outDir=None, classifier=None, classifierArgs=None, develFolds=5, useCafa=False, useCombinations=True, useLearning=True, baselineCutoff=1, clear=False):
    if outDir != None:
        if clear and os.path.exists(outDir):
            print "Removing output directory", outDir
            shutil.rmtree(outDir)
        if not os.path.exists(outDir):
            print "Making output directory", outDir
            os.makedirs(outDir)
        Stream.openLog(os.path.join(outDir, "log.txt"))
    
    print "==========", "Ensemble", "=========="
    proteins = {}
    print "Loading Swissprot proteins"
    loading.loadFASTA(os.path.join(dataPath, "Swiss_Prot", "Swissprot_sequence.tsv.gz"), proteins)
    if useCafa == "cafa":
        print "Loading CAFA3 targets"
        loading.loadFASTA(os.path.join(dataPath, "CAFA3_targets", "Target_files", "target.all.fasta"), proteins, True)
    print "Proteins:", len(proteins)
    termCounts = loading.loadAnnotations(os.path.join(dataPath, "data", "Swissprot_propagated.tsv.gz"), proteins)
    print "Unique terms:", len(termCounts)
    topTerms = loading.getTopTerms(termCounts, 5000)
    limitTerms=set([x[0] for x in topTerms])
    print "Using", len(topTerms), "most common GO terms"
    loading.loadSplit(os.path.join(dataPath, "data"), proteins)
    loading.defineSets(proteins, "overlap" if useCafa else "skip")
    
    if nnInput != None:
        print "Loading neural network predictions from", nnInput
        for setName in (("devel", "test", "cafa") if useCafa else ("devel", "test")):
            evaluateFile.loadPredictions(proteins, os.path.join(nnInput, setName + ("_targets" if setName == "cafa" else "_pred")) + ".tsv.gz", limitToSets=None, readGold=False, predKey="nn_pred", confKey="nn_pred_conf")
    if clsInput != None:
        print "Loading classifier predictions"
        evaluateFile.loadPredictions(proteins, clsInput, limitToSets=["devel","test","cafa"] if useCafa else ["devel","test"], readGold=True, predKey="cls_pred", confKey="cls_pred_conf")
    if baselineCutoff > 0:
        print "Loading baseline predictions"
        loading.loadBaseline(dataPath, proteins, "baseline_pred", baselineCutoff, topTerms)
    
    if useCombinations:
        print "===============", "Combining predictions", "===============" 
        combKey = "combined"
        predKeys = ["nn_pred", "cls_pred"]
        if nnInput != None:
            predKeys += ["nn_pred"]
        if clsInput != None:
            predKeys += ["cls_pred"]
        if baselineCutoff > 0:
            predKeys += ["baseline_pred"]
        combinations = getCombinations(predKeys)
        print "Testing", len(combinations), "combinations"
        for combination in combinations:
            print "******", "Combination", combination, "******"
            for setName in ("devel", "test"):
                for mode in (("AND", "OR") if len(combination) > 1 else ("SINGLE",)):
                    print "***", "Evaluating predictions for set '" + setName + "' using mode '" + mode + "'", "***"
                    combinePred(proteins, combination, combKey, mode, limitToSets=[setName])
                    examples = evaluateFile.makeExamples(proteins, limitTerms=limitTerms, limitToSets=[setName], predKey=combKey)
                    loading.vectorizeExamples(examples, None)
                    results = evaluation.evaluate(examples["labels"], examples["predictions"], examples, terms=None, averageOnly=True)
                    print "Average for", str(combination) + "/" + setName + "/" + mode + ":", evaluation.metricsToString(results["average"])
    
    if useLearning:
        print "===============", "Learning", "==============="
        Classifier = classification.importNamed(classifier)
        examples = buildExamples(proteins, "nn_pred", "cls_pred", limitToSets=["devel", "test"], limitTerms=limitTerms, outDir=outDir)
        learn(examples, Classifier, classifierArgs, develFolds=develFolds, limitTerms=limitTerms)

if __name__=="__main__":       
    from optparse import OptionParser
    optparser = OptionParser(description="Combine relation predictions (All input files must include both positive and negative interaction elements)")
    optparser.add_option("-p", "--dataPath", default=os.path.expanduser("~/data/CAFA3/data"), help="")
    optparser.add_option("-a", "--nnInput", default=None)
    optparser.add_option("-b", "--clsInput", default=None)
    optparser.add_option("-g", "--gold", default=None)
    optparser.add_option("-o", "--outDir", default=None)
    optparser.add_option("-s", "--simple", default=False, action="store_true")
    optparser.add_option("-l", "--learning", default=False, action="store_true")
    optparser.add_option("-f", "--baseline", default=-1, type=int)
    #optparser.add_option("-w", "--write", default="OR")
    optparser.add_option("-n", "--develFolds", type=int, default=5)
    optparser.add_option('-c','--classifier', help='', default="ensemble.RandomForestClassifier")
    optparser.add_option('-r','--args', help='', default="{'random_state':[1], 'n_estimators':[10], 'n_jobs':[1], 'verbose':[3]}")
    optparser.add_option("--clear", default=False, action="store_true", help="Remove the output directory if it already exists")
    (options, args) = optparser.parse_args()
    
    #assert options.write in ("AUTO", "LEARN", "AND", "OR")
    options.args = eval(options.args)
    combine(dataPath=options.dataPath, nnInput=options.nnInput, clsInput=options.clsInput, outDir=options.outDir,
            classifier=options.classifier, classifierArgs=options.args, develFolds=options.develFolds,
            useCombinations=options.simple, useLearning=options.learning, baselineCutoff=options.baseline,
            clear=options.clear)