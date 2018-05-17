import sys, os
import learning.loading as loading
import learning.makeFolds as makeFolds
import operator
from collections import Counter
from learning.featureBuilders import *

class Task(object):
    def __init__(self, dataPath):
        self.proteins = None
        
        self.dataPath = dataPath
        self.sequencesPath = None
        self.targetsPath = None
        self.annotationsPath = None
        self.splitPath = None
        self.foldsPath = None
        
        self.builders = {"taxonomy":TaxonomyFeatureBuilder,
                         }
        
        self.features = {"taxonomy":TaxonomyFeatureBuilder(["Taxonomy"]),
                         "similar":UniprotFeatureBuilder(os.path.join(dataPath, "Uniprot", "similar.txt")),}
        
        self.numTerms = 1000
        self.removeNonHuman = False
        self.remapSets = None
        self.allowMissing = False
        self.limitTrainingToAnnotated = False
    
    def getPath(self, subPathComponents):
        return os.path.join(self.dataPath, subPathComponents)
    
    def loadProteins(self, cafaTargets="skip"):
        assert cafaTargets in ("skip", "overlap", "separate", "external")
        self.cafaTargets = cafaTargets
        self.proteins = {}
        loading.loadFASTA(self.getPath(self.sequencesPath), self.proteins)
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
    
    def buildFeatures(self):
        featuresDataPath = dataPath
        if task == "cafapi":
            featuresDataPath = os.path.join(dataPath, "CAFA_PI", "features")
        print "Loading features from", featuresDataPath
        examples = buildExamples(proteins, featuresDataPath, limit, limitTerms=set([x[0] for x in topTerms]), featureGroups=featureGroups, debug=debug)
        print "Saving examples to", exampleFilePath
        with gzip.open(exampleFilePath, "wt") as pickleFile:
            json.dump(examples, pickleFile, indent=2) #pickle.dump(examples, pickleFile)
        loading.vectorizeExamples(examples, featureGroups)
    
    def getTopTerms(self, counts, num=1000):
        return sorted(counts.items(), key=operator.itemgetter(1), reverse=True)[0:num]
        

class CAFAPITask(Task):
    def __init__(self, dataPath):
        Task.__init__(self, dataPath)
        
        self.sequencesPath = self.getPath(["CAFA_PI", "Swissprot", "CAFA_PI_Swissprot_sequence.tsv.gz"])
        self.targetsPath = self.getPath(["CAFA_PI", "Swissprot", "target.all.fasta.gz"])
        self.annotationsPath = self.getPath(["CAFA_PI", "Swissprot", "CAFA_PI_Swissprot_propagated.tsv.gz"])
        self.splitPath = self.getPath(["CAFA_PI", "Swissprot"])
        self.foldsPath = self.getPath(["folds", "CAFA_PI_training_folds_180417.tsv.gz"])
        
        self.features = {"taxonomy":TaxonomyFeatureBuilder(["Taxonomy"]),
                         "similar":UniprotFeatureBuilder(os.path.join(dataPath, "Uniprot", "similar.txt")),}