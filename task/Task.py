import sys, os
import gzip
import learning.loading as loading
import learning.makeFolds as makeFolds
import operator
from collections import Counter
import json
from learning.classification import Classification, SingleLabelClassification
import utils.statistics as statistics

class Task(object):
    ###########################################################################
    # Task Definitions
    ###########################################################################
    TASKS = {}
    
    @staticmethod
    def registerTask(task, name):
        assert name not in Task.TASKS
        Task.TASKS[name] = task
    
    @staticmethod
    def getTask(name):
        assert name in Task.TASKS
        return Task.TASKS[name]()
    
    ###########################################################################
    # Initialization
    ###########################################################################
    
    def __init__(self):
        # Internal data structures
        self.proteins = None
        self.examples = None
        self.debug = False
        
        # Task Specific Settings ##############################################
        
        # File paths
        self.dataPath = None # The main data path (all other paths are relative to this)
        self.sequencesPath = None # Sequences in the FASTA format
        self.targetsPath = None # Sequences in the FASTA format
        self.annotationsPath = None # Ontology term annotations for the sequences
        self.splitPath = None # The directory containing the train/devel/test split
        self.foldsPath = None # A file containing the n-fold cross-validation groups
        # Feature Groups
        self.features = None # A dictionary of feature group name / FeatureBuilder pairs
        self.defaultFeatures = None # The list of the names of the feature groups which are used by default
        # Task Configuration
        self.numTerms = 5000 # The cutoff for the n most common terms
        self.removeNonHuman = False # Remove all proteins where the species != human (for the HPO task)
        self.remapSets = None # A dictionary for remapping sets, e.g. {"test":"devel"}
        self.allowMissing = False # Whether all protein ids in the train/devel/set must exists among the loaded proteins
        self.limitTrainingToAnnotated = False # Remove all proteins which have no annotated terms
        self.annotationFormat = "GO" # 'GO' or 'HPO' annotation file format
    
    ###########################################################################
    # Loading
    ###########################################################################
    
    def setDebug(self, debug):
        self.debug = debug
        for group in self.features:
            self.features[group].debug = debug
    
    def setDataPath(self, dataPath):
        assert self.dataPath == None
        self.dataPath = dataPath
        if dataPath != None:
            self.sequencesPath = self._getPath(self.sequencesPath)
            self.targetsPath = self._getPath(self.targetsPath)
            self.annotationsPath = self._getPath(self.annotationsPath)
            self.splitPath = self._getPath(self.splitPath)
            self.foldsPath = self._getPath(self.foldsPath)
            self.termsPath = self._getPath(self.termsPath)
            for group in self.features:
                self.features[group].setDataPath(dataPath)
    
    def _getPath(self, subPath):
        return os.path.join(self.dataPath, subPath) if subPath != None else None
    
    def loadProteins(self, cafaTargets="skip"):
        assert cafaTargets in ("skip", "overlap", "separate", "external")
        self.cafaTargets = cafaTargets
        self.proteins = {}
        loading.loadFASTA(self.sequencesPath, self.proteins)
        if cafaTargets != "skip" and self.targetsPath != None:
            loading.loadFASTA(self.targetsPath, self.proteins, True)
        if self.removeNonHuman:
            loading.removeNonHuman(self.proteins)
        assert self.annotationFormat in ("GO", "HPO")
        if self.annotationsPath:
            if self.annotationFormat == "GO":
                self.termCounts = loading.loadAnnotations(self.annotationsPath, self.proteins)
            else:
                self.termCounts = loading.loadHPOAnnotations(self.annotationsPath, self.proteins)
        else:
            print "No annotations to load"
            self.termCounts = {}
        print "Unique terms:", len(self.termCounts)
        topTerms = self.getTopTerms(self.termCounts, self.numTerms)
        print "Using", len(topTerms), "most common GO terms"
    
    def getTopTerms(self, counts, num=1000):
        return sorted(counts.items(), key=operator.itemgetter(1), reverse=True)[0:num]
    
    def loadSplit(self, fold=None):
        if self.splitPath == None:
            print "No split to load"
            for protId in self.proteins.keys():
                self.proteins[protId]["split"] = "undefined"
        else:
            loading.loadSplit(self.splitPath, self.proteins, self.allowMissing)

        if self.remapSets != None:
            for protId in self.proteins.keys():
                protSet = self.proteins[protId].get("split")
                if protSet in self.remapSets:
                    self.proteins[protId]["split"] = self.remapSets[protSet]
            print "Remapped splits", Counter([x.get("split") for x in self.proteins.values()])
        if fold != None:
            makeFolds.loadFolds(self.proteins, self.foldsPath)
        
        loading.defineSets(self.proteins, self.cafaTargets, fold=fold, limitTrainingToAnnotated = self.limitTrainingToAnnotated)
    
    ###########################################################################
    # Example Generation
    ###########################################################################
    
    def _getFeatureGroups(self, groups=None):
        if groups == None:
            groups = self.defaultFeatures if self.defaultFeatures != None else sorted(self.features.keys())
        if "all" in groups:
            groups.remove("all")
            groups = groups + sorted(self.features.keys())
        removed = [x for x in groups if x.startswith("-")]
        groups = [x for x in groups if not x.startswith("-")]
        for group in removed:
            groups.remove(group.strip("-"))
        return groups
    
    def buildExamples(self, groups=None, limit=None, limitTerms=None, featureGroups=None):
        print "Building examples"
        self.examples = {"labels":[], "features":[], "ids":[], "cafa_ids":[], "sets":[], "label_names":[], "label_size":{}}
        protIds = sorted(self.proteins.keys())
        if limit:
            protIds = protIds[0:limit]
        protObjs = [self.proteins[key] for key in protIds]
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
                if label not in self.examples["label_size"]:
                    self.examples["label_size"][label] = 0
                self.examples["label_size"][label] += 1
            self.examples["labels"].append(labels)
            self.examples["ids"].append(protein["id"])
            self.examples["cafa_ids"].append(protein["cafa_ids"])
            self.examples["sets"].append(protein["sets"])
        # Build features
        groups = self._getFeatureGroups(groups)
        print "Building features, feature groups =", featureGroups
        for group in groups:
            if group not in self.features.keys():
                raise Exception("Unknown feature group '" + str(group) + "'")
            print "Building features for group", group
            self.features[group].build(protObjs)
        self.examples["features"] = [x["features"] for x in protObjs]
        for protObj in protObjs:
            del protObj["features"]
        # Prepare the examples
        print "Built", len(self.examples["labels"]), "examples" # with", len(examples["feature_names"]), "unique features"
        return self.examples
    
    def saveExamples(self, exampleFilePath):
        print "Saving examples to", exampleFilePath
        with gzip.open(exampleFilePath, "wt") as pickleFile:
            json.dump(self.examples, pickleFile, indent=2, sort_keys=True)
    
    def loadExamples(self, exampleFilePath):
        print "Loading examples from", exampleFilePath
        with gzip.open(exampleFilePath, "rt") as pickleFile:
            self.examples = json.load(pickleFile)
    
    ###########################################################################
    # Classification
    ###########################################################################
    
    def vectorizeExamples(self, idPath=None):
        loading.vectorizeExamples(self.examples, idPath=idPath)
    
    def train(self, outDir, classifier=None, classifierArgs=None, singleLabelJobs=None, negatives=False, useTestSet=False):
        terms = loading.loadGOTerms(self.termsPath)
        loading.saveIdNames(self.examples["feature_names"], os.path.join(outDir, "features.tsv"))
        loading.saveIdNames(self.examples["label_names"], os.path.join(outDir, "labels.tsv"))
        if singleLabelJobs == None:
            cls = Classification()
        else:
            cls = SingleLabelClassification(singleLabelJobs)
        cls.optimize(classifier, classifierArgs, self.examples, terms=terms, 
                     outDir=outDir, negatives=negatives,
                     useTestSet=useTestSet, useCAFASet=(self.cafaTargets != "skip"))
    
    def classify(self, outDir, modelPath, singleLabelJobs=None, negatives=False, useTestSet=False):
        terms = loading.loadGOTerms(self.termsPath)
        loading.saveIdNames(self.examples["feature_names"], os.path.join(outDir, "features.tsv"))
        loading.saveIdNames(self.examples["label_names"], os.path.join(outDir, "labels.tsv"))
        if singleLabelJobs == None:
            cls = Classification()
        else:
            cls = SingleLabelClassification(singleLabelJobs)
        cls.predict(modelPath, self.examples, terms=terms, outDir=outDir, negatives=negatives, useTestSet=useTestSet, useCAFASet=(self.cafaTargets != "skip"))
    
    def makeStatistics(self, outDir):
        statistics.makeStatistics(self.examples, outDir)