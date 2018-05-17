import sys, os
import learning.loading as loading
import learning.makeFolds as makeFolds
import operator
from collections import Counter
from learning.featureBuilders import *
import json

class Task(object):
    def __init__(self):
        self.proteins = None
        self.examples = None
        
        self.dataPath = None
        self.sequencesPath = None
        self.targetsPath = None
        self.annotationsPath = None
        self.splitPath = None
        self.foldsPath = None
        
        self.features = None
        self.defaultFeatures = None
        
        self.numTerms = 1000
        self.removeNonHuman = False
        self.remapSets = None
        self.allowMissing = False
        self.limitTrainingToAnnotated = False
    
    def setDataPath(self, dataPath):
        assert self.dataPath == None
        self.dataPath = dataPath
        if dataPath != None:
            self.sequencesPath = self._getPath(self.sequencesPath)
            self.targetsPath = self._getPath(self.targetsPath)
            self.annotationsPath = self._getPath(self.annotationsPath)
            self.splitPath = self._getPath(self.splitPath)
            self.foldsPath = self._getPath(self.foldsPath)
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
        self.termCounts = loading.loadAnnotations(self.annotationsPath, self.proteins)
        print "Unique terms:", len(self.termCounts)
        topTerms = self.getTopTerms(self.termCounts, self.numTerms)
        print "Using", len(topTerms), "most common GO terms"
    
    def loadSplit(self, fold=None):
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
    
    def _getFeatureGroups(self, groups=None):
        if groups == None:
            groups = self.defaultFeatures if self.defaultFeatures != None else sorted(self.features.keys())
            groups = ["all"]
        return [x for x in groups if not x.startswith("-")]
    
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
            json.dump(self.examples, pickleFile, indent=2)
    
    def loadExamples(self, exampleFilePath):
        print "Loading examples from", exampleFilePath
        with gzip.open(exampleFilePath, "rt") as pickleFile:
            self.examples = json.load(pickleFile)
    
    def vectorizeExamples(self):
        loading.vectorizeExamples(self.examples)
    
    def getTopTerms(self, counts, num=1000):
        return sorted(counts.items(), key=operator.itemgetter(1), reverse=True)[0:num]
        

class CAFAPITask(Task):
    def __init__(self):
        Task.__init__(self)
        
        self.sequencesPath = "CAFA_PI/Swissprot/CAFA_PI_Swissprot_sequence.tsv.gz"
        self.targetsPath = "CAFA_PI/Swissprot/target.all.fasta.gz"
        self.annotationsPath = "CAFA_PI/Swissprot/CAFA_PI_Swissprot_propagated.tsv.gz"
        self.splitPath = "CAFA_PI/Swissprot"
        self.foldsPath = "folds/CAFA_PI_training_folds_180417.tsv.gz"
        
        self.features = {
            "taxonomy":TaxonomyFeatureBuilder(["Taxonomy"]),
            "similar":UniprotFeatureBuilder("Uniprot/similar.txt"),
            "blast":BlastFeatureBuilder(["temp_blastp_result_features", "blastp_result_features"]),
            "blast62":BlastFeatureBuilder(["CAFA2/training_features", "CAFA2/CAFA3_features"], tag="BLAST62"),
            "delta":BlastFeatureBuilder(["temp_deltablast_result_features", "deltablast_result_features"], tag="DELTA"),
            "interpro":InterProScanFeatureBuilder(["temp_interproscan_result_features", "interproscan_result_features"]),
            "predgpi":PredGPIFeatureBuilder(["predGPI"]),
            "nucpred":NucPredFeatureBuilder(["nucPred"]),
            "netacet":NetAcetFeatureBuilder(["NetAcet"]),
            "funtaxis":FunTaxISFeatureBuilder(["FunTaxIS"]),
            "ngrams":NGramFeatureBuilder(["ngrams/4jari/min_len3-min_freq2-min1fun-top_fun5k"])
        }