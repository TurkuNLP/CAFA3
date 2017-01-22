from evaluation import evaluate, metricsToString, getResultsString, resultIsBetter, saveResults
from sklearn.grid_search import ParameterGrid, GridSearchCV
import os
from sklearn.preprocessing.label import LabelBinarizer
from sklearn.externals.joblib.parallel import Parallel, delayed
from sklearn.multiclass import _fit_binary

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
        
    def getSubset(self, examples, setNames):
        sets = examples["sets"]
        indices = [i for i in range(len(sets)) if any(x in setNames for x in sets[i])]
        features = examples["features"][indices]
        labels = examples["labels"][indices]
        ids = [examples["ids"][i] for i in range(len(sets)) if any(x in setNames for x in sets[i])]
        cafa_ids = [examples["cafa_ids"][i] for i in range(len(sets)) if any(x in setNames for x in sets[i])]
        return features, labels, indices, ids, cafa_ids
    
    def learn(self, args, examples, trainSets, testSets, terms, cls=None, averageOnly=False):
        args = args.copy()
        print "Learning with args", args
        if cls == None:
            print "Initializing classifier", self.Classifier.__name__, "with arguments", args
            cls = self.Classifier(**args)
        else:
            print "Using existing classifier"
        print "Train / test = ", trainSets, "/", testSets
        trainFeatures, trainLabels, _, _, _ = self.getSubset(examples, trainSets)
        testFeatures, testLabels, _, testIds, testCafaIds = self.getSubset(examples, testSets)
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
            cls, data = self.learn(args, examples, ["train"], ["devel"], False, terms, cls=cls, averageOnly=True)
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
    
    def predictSet(self, args, examples, trainSets, testSets, terms, outDir):
        print "Classifying set", testSets, "using sets", trainSets
        _, data = self.learn(best["args"], examples, ["train", "devel"], ["test"], terms)
        if outDir != None:
            saveResults(data, os.path.join(outDir, "test"), examples["label_names"], negatives=negatives)
    
    def optimize(self, classifier, classifierArgs, examples, terms=None, outDir=None, negatives=False, useTestSet=False, useCAFASet=False):
        best = None
        print "Parameter grid search"
        self.Classifier = importNamed(classifier)
        if classifierArgs.get("warm_start") == [True]:
            best = self.warmStartGrid(classifierArgs, examples, terms)
        else:
            for args in ParameterGrid(classifierArgs):
                _, data = self.learn(args, examples, ["train"], ["devel"], terms)
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
            _, data = self.learn(best["args"], examples, ["train", "devel"], ["test"], terms)
            if outDir != None:
                saveResults(data, os.path.join(outDir, "test"), examples["label_names"], negatives=negatives)
        if useCAFASet:
            print "Classifying the CAFA targets"
            _, data = self.learn(best["args"], examples, ["train", "devel", "test"], ["cafa"], terms)
            if outDir != None:
                saveResults(data, os.path.join(outDir, "cafa"), examples["label_names"], negatives=negatives)
                
class SingleLabelClassification(Classification):
    def __init__(self, n_jobs):
        self.n_jobs = n_jobs
    
    def optimize(self, classifier, classifierArgs, examples, terms=None, outDir=None, negatives=False, useTestSet=False, useCAFASet=False):
        print "Parameter grid search"
        self.Classifier = importNamed(classifier)
        best = None
        examples["multilabels"] = examples["labels"]
        for labelIndex in range(len(examples["label_names"])):
            print "Parameter search for label", labelIndex, examples["label_names"][labelIndex]
            examples["labels"] = examples["multilabels"][:, labelIndex]
            trainFeatures, trainLabels, _, _, _ = self.getSubset(examples, ["train"])
            testFeatures, testLabels, _, testIds, testCafaIds = self.getSubset(examples, ["devel"])
            clf = GridSearchCV(self.Classifier, classifierArgs, "f1score", n_jobs=self.n_jobs)
            clf.fit(trainFeatures, trainLabels)
            print "Best params", (clf.best_params_, clf.best_score_)
            print "Predicting"
            testPredictions = clf.predict(testFeatures)
        
        if classifierArgs.get("warm_start") == [True]:
            best = self.warmStartGrid(classifierArgs, examples, terms)
        else:
            for args in ParameterGrid(classifierArgs):
                _, data = self.learn(args, examples, ["train"], ["devel"], terms)
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
            _, data = self.learn(best["args"], examples, ["train", "devel"], ["test"], terms)
            if outDir != None:
                saveResults(data, os.path.join(outDir, "test"), examples["label_names"], negatives=negatives)
        if useCAFASet:
            print "Classifying the CAFA targets"
            _, data = self.learn(best["args"], examples, ["train", "devel", "test"], ["cafa"], terms)
            if outDir != None:
                saveResults(data, os.path.join(outDir, "cafa"), examples["label_names"], negatives=negatives)