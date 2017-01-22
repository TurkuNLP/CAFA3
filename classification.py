from evaluation import evaluate, metricsToString, getResultsString

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

class Classification():
    def __init__(self):
        self.Classifier = None
    
    def learn(self, args, examples, trainSets, testSets, terms, cls=None, averageOnly=False):
        args = args.copy()
        print "Learning with args", args
        if cls == None:
            print "Initializing classifier", self.Classifier.__name__, "with arguments", args
            cls = self.Classifier(**args)
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
    
    def resultIsBetter(self, original, new, key="average"):
        if new[key]["fscore"] != original[key]["fscore"]:
            return new[key]["fscore"] > original[key]["fscore"]
        else:
            return new[key]["auc"] > original[key]["auc"]
    
    def warmStartGrid(self, classifierArgs, examples, terms):
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
    
    def optimize(self, classifier, classifierArgs, examples, cvJobs=1, terms=None, useMultiOutputClassifier=False, outDir=None, negatives=False, useTestSet=False, useCAFASet=False):
        best = None
        print "Parameter grid search"
        self.Classifier = importNamed(classifier)
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