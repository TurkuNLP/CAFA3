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

def clearKeys(proteins, keys):
    for protId in proteins:
        protein = proteins[protId]
        for key in keys:
            if key in protein:
                del protein[key]

def getCombinations(items):
    combinations = []
    for i in xrange(1, len(items) + 1):
        els = [list(x) for x in itertools.combinations(items, i)]
        combinations.extend(els)
    return combinations

def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)

def combineConf(protein, labels, predKeys, combKey):
    combConfKey = combKey + "_conf"
    assert combConfKey not in protein
    protein[combConfKey] = {}
    combConfs = protein[combConfKey]
    assert combKey + "_sources" not in protein
    protein[combKey + "_sources"] = {}
    combSources = protein[combKey + "_sources"]
    for key in predKeys:
        preds = protein.get(key, {})
        predConfs = protein.get(key + "_conf", {})
        for label in labels:
            predConf = predConfs.get(label)
            if predConf != None:
                if label not in combConfs:
                    combConfs[label] = []
                combConfs[label].append(predConf)
            if label in preds:
                if label not in combSources:
                    combSources[label] = []
                combSources[label].append(key)
    protein[combConfKey] = {x:mean(combConfs[x]) for x in combConfs}
                        
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
        assert combKey not in protein
        protein[combKey] = {}
        if mode == "OR" or mode == "SINGLE":
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
                protein[combKey] = {x:1 for x in set.intersection(*predLabelSets)} #{x:1 for x in pred1.intersection(pred2)}
        combineConf(protein, sorted(protein[combKey].keys()), predKeys, combKey)            
    print "Combined predictions, mode =", mode, "counts =", dict(counts)

def buildFeatures(protein, label, predKeys, predictions, confidences, counts):
    features = {label:1}
    for key in predKeys:
        if key in predictions:
            features["pos:" + key] = 1
            if key in confidences:
                if label in confidences[key]:
                    features["conf:" + key] = confidences[key].get(label)
        else:
            features["neg:" + key] = 1
    return features

def buildExamples(proteins, predKeys, limitToSets=None, limitTerms=None, outDir=None):
    counts = defaultdict(int)
    empty = {}
    examples = {"classes":[], "features":[], "sets":[], "proteins":[], "labels":[]}
    for protId in proteins:
        protein = proteins[protId]
        if limitToSets != None and not any(x in limitToSets for x in protein["sets"]):
            counts["out-of-sets"] += 1
            continue
        predictions = {}
        confidences = {}
        predLabels = []
        for key in predKeys: # Collect predictions for each system (key)
            if key not in protein:
                counts["no-prediction-for-" + key] += 1
            else:
                counts["predictions-for-" + key] += 1
                predLabels.extend(protein[key].keys())
                predictions[key] = protein[key]
                confidences[key] = protein.get(key + "_conf", empty)
        predLabels = sorted(set(predLabels))
        #predLabels = sorted(set.union(*[predictions[x].keys() for x in predictions])) # + goldLabels))
        protSets = protein.get("sets")
        goldLabels = protein["terms"].keys()
        if limitTerms:
            predLabels = [x for x in predLabels if x in limitTerms]
            goldLabels = [x for x in goldLabels if x in limitTerms]
        goldLabels = set(goldLabels)
        for label in predLabels:
            features = buildFeatures(protein, label, predKeys, predictions, confidences, counts)
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

def binaryToMultiLabel(examples, predictions, probabilities, predKey):
    print "Converting binary predictions to labels"
    predKeyConf = predKey + "_conf"
    for prediction, probability, protein, label in zip(predictions, probabilities, examples["proteins"], examples["labels"]):
        if predKey not in protein:
            protein[predKey] = {}
            protein[predKeyConf] = {}
        if prediction == 1:
            protein[predKey][label] = 1
            protein[predKeyConf][label] = probability

def learn(examples, Classifier, classifierArgs, develFolds=10, verbose=3, n_jobs=1, predKey="ml_comb_pred", limitTerms=None):
    print "Parameter grid search"
    develExamples = getSubset(examples, ["devel"])
    clf = GridSearchCV(Classifier(), classifierArgs, cv=develFolds, verbose=verbose, n_jobs=n_jobs, scoring="f1_micro")
    clf.fit(develExamples["features"], develExamples["classes"])
    print "Best params", (clf.best_params_, clf.best_score_)
    print "Predicting the test set"
    testExamples = getSubset(examples, ["test"])
    testPredictions = clf.predict(testExamples["features"])
    testProbabilities = clf.predict_proba(testExamples["features"])
    binaryToMultiLabel(testExamples, testPredictions, testProbabilities, predKey)
    print "Evaluating test set ensemble predictions"
    testProteins = {x["id"]:x for x in testExamples["proteins"]}
    multiLabelTestExamples = evaluateFile.makeExamples(testProteins, limitTerms=limitTerms, limitToSets=["test"], predKey=predKey)
    loading.vectorizeExamples(multiLabelTestExamples, None, sparseLabels=True)
    results = evaluation.evaluate(multiLabelTestExamples["labels"], multiLabelTestExamples["predictions"], multiLabelTestExamples, terms=None, averageOnly=True, noAUC=True)
    print "Average for test set:", evaluation.metricsToString(results["average"])
    print "Predicting all examples"
    allPredictions = clf.predict(examples["features"])
    allProbabilities = clf.predict_proba(examples["features"])
    binaryToMultiLabel(testExamples, allPredictions, allProbabilities, predKey)
    
def combine(dataPath, nnInput, clsInput, outDir=None, classifier=None, classifierArgs=None, develFolds=5, useCafa=False, useCombinations=True, useLearning=True, baselineCutoff=1, numTerms=5000, clear=False, useOutFiles=True):
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
    if useCafa:
        print "Loading CAFA3 targets"
        loading.loadFASTA(os.path.join(dataPath, "CAFA3_targets", "Target_files", "target.all.fasta"), proteins, True)
    print "Proteins:", len(proteins)
    termCounts = loading.loadAnnotations(os.path.join(dataPath, "data", "Swissprot_propagated.tsv.gz"), proteins)
    print "Unique terms:", len(termCounts)
    topTerms = loading.getTopTerms(termCounts, numTerms)
    limitTerms=set([x[0] for x in topTerms])
    print "Using", len(topTerms), "most common GO terms"
    loading.loadSplit(os.path.join(dataPath, "data"), proteins)
    loading.defineSets(proteins, "overlap" if useCafa else "skip")
    
    predKeys = []
    if nnInput != None:
        print "Loading neural network predictions from", nnInput
        for setName in (("devel", "test", "cafa") if useCafa else ("devel", "test")):
            predKey = "nn_pred_cafa" if setName == "cafa" else "nn_pred"
            evaluateFile.loadPredictions(proteins, os.path.join(nnInput, setName + ("_targets" if setName == "cafa" else "_pred")) + ".tsv.gz", limitToSets=None, readGold=False, predKey=predKey, confKey="nn_pred_conf", includeDuplicates=True)
        predKeys += ["nn_pred"]
    if clsInput != None:
        print "Loading classifier predictions"
        evaluateFile.loadPredictions(proteins, clsInput, limitToSets=["devel","test","cafa"] if useCafa else ["devel","test"], readGold=True, predKey="cls_pred", confKey="cls_pred_conf")
        predKeys += ["cls_pred"]
    if baselineCutoff > 0:
        print "Loading baseline predictions"
        loading.loadBaseline(dataPath, proteins, "baseline_pred", baselineCutoff, limitTerms, useCafa=useCafa)
        predKeys += ["baseline_pred"]
    
    if useCombinations:
        print "===============", "Combining predictions", "===============" 
        combKey = "comb_pred"
        combConfKey = "comb_pred_conf"
        combinations = getCombinations(predKeys)
        numCombinations = len(combinations)
        print "Testing", numCombinations, "combinations"
        for i in range(len(combinations)):
            print "******", "Combination", str(i + 1) + "/" + str(numCombinations), combinations[i], "******"
            for mode in (("AND", "OR") if len(combinations[i]) > 1 else ("SINGLE",)):
                for setName in (("devel", "test", "cafa") if useCafa else ("devel", "test")):
                    combination = combinations[i][:]
                    if setName == "cafa" and "nn_pred" in combination:
                        combination[combination.index("nn_pred")] = "nn_pred_cafa"
                    print "***", "Evaluating", combination, "predictions for set '" + setName + "' using mode '" + mode + "'", "***"
                    combinePred(proteins, combination, combKey, mode, limitToSets=[setName])
                    if setName != "cafa":
                        examples = evaluateFile.makeExamples(proteins, limitTerms=limitTerms, limitToSets=[setName], predKey=combKey)
                        loading.vectorizeExamples(examples, None, sparseLabels=True)
                        results = evaluation.evaluate(examples["labels"], examples["predictions"], examples, terms=None, averageOnly=True, noAUC=True)
                        print "Average for", str(combination) + "/" + setName + "/" + mode + ":", evaluation.metricsToString(results["average"])
                    else:
                        print "Skipping evaluation for set '" + setName + "'"
                    if useOutFiles:
                        combString = "-".join(combination)
                        outPath = os.path.join(outDir, "-".join([combString, setName, mode, "ensemble"]) + ".tsv.gz")
                        evaluation.saveProteins(proteins, outPath, limitTerms=limitTerms, limitToSets=[setName], predKey=combKey) #pass#evaluation.saveResults(data, outStem, label_names, negatives)
                    clearKeys(proteins, [combKey, combConfKey, combKey + "_sources"])
    if useLearning:
        print "===============", "Learning", "==============="
        Classifier = classification.importNamed(classifier)
        examples = buildExamples(proteins, predKeys, limitToSets=None, limitTerms=limitTerms, outDir=outDir)
        learn(examples, Classifier, classifierArgs, develFolds=develFolds, limitTerms=limitTerms, predKey="ml_comb_pred")
        if useOutFiles:
            for setName in (("devel", "test", "cafa") if useCafa else ("devel", "test")):
                outPath = os.path.join(outDir, "-".join([setName, "ML", "ensemble"]) + ".tsv.gz")
                evaluation.saveProteins(proteins, outPath, limitTerms=limitTerms, limitToSets=[setName], predKey="ml_comb_pred") #pass#evaluation.saveResults(data, outStem, label_names, negatives)

if __name__=="__main__":       
    from optparse import OptionParser
    optparser = OptionParser(description="Ensemble")
    optparser.add_option("-p", "--dataPath", default=os.path.expanduser("~/data/CAFA3/data"), help="Data directory")
    optparser.add_option("-a", "--nnInput", default=None, help="Neural network predictions tsv.gz file")
    optparser.add_option("-b", "--clsInput", default=None, help="Classifier predictions tsv.gz file")
    optparser.add_option("-o", "--outDir", default=None, help="Output directory")
    optparser.add_option("-s", "--simple", default=False, action="store_true", help="Do simple ensembles (AND, OR and SINGLE)")
    optparser.add_option("-l", "--learning", default=False, action="store_true", help="Do machine learning ensemble")
    optparser.add_option("-f", "--baseline", default=-1, type=int, help="Add the BLAST baseline as the third set of predictions. Value in range 1-10, 1 for all values.")
    optparser.add_option("-t", "--terms", default=5000, type=int, help="The number of top most common GO terms to use as labels")
    optparser.add_option("-w", "--write", default=False, action="store_true", help="Write output files")
    optparser.add_option("-n", "--develFolds", type=int, default=5, help="Cross-validation for parameter optimization")
    optparser.add_option('-c','--classifier', default="ensemble.RandomForestClassifier", help="Scikit-learn classifier")
    optparser.add_option('-r','--args', default="{'random_state':[1], 'n_estimators':[10], 'n_jobs':[1], 'verbose':[3]}", help="Classifier arguments")
    optparser.add_option("--clear", default=False, action="store_true", help="Remove the output directory if it already exists")
    optparser.add_option("--cafa", default=False, action="store_true", help="Process CAFA predictions")
    (options, args) = optparser.parse_args()
    
    options.args = eval(options.args)
    combine(dataPath=options.dataPath, nnInput=options.nnInput, clsInput=options.clsInput, outDir=options.outDir,
            classifier=options.classifier, classifierArgs=options.args, develFolds=options.develFolds,
            useCafa=options.cafa,
            useCombinations=options.simple, useLearning=options.learning, baselineCutoff=options.baseline,
            numTerms=options.terms, clear=options.clear, useOutFiles=options.write)