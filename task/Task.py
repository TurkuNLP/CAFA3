import sys, os
import learning.loading as loading
import operator
from collections import Counter

class Task:
    def __init__(self):
        self.proteins = {}
        self.settings = {"numTerms":1000}
        self.removeNonHuman = False
        self.remapSets = None
    
    def load(self, cafaTargets="skip"):
        assert cafaTargets in ("skip", "overlap", "separate", "external")
        self.proteins = {}
        loading.loadFASTA(self.settings["sequences"], self.proteins)
        if cafaTargets != "skip" and self.settings.get("targets") != None:
            loading.loadFASTA(self.settings["targets"], self.proteins, True)
        if self.removeNonHuman:
            loading.removeNonHuman(self.proteins)
        self.termCounts = loading.loadAnnotations(self.settings.get("annotations"), self.proteins)
        print "Unique terms:", len(self.termCounts)
        topTerms = self.getTopTerms(self.termCounts, self.settings.get("numTerms"))
        print "Using", len(topTerms), "most common GO terms"
    
    def loadSplit(self):
        loading.loadSplit(self.settings["split"], self.proteins, self.allowMissing)
        if self.remapSets != None:
            for protId in self.proteins.keys():
                protSet = self.proteins[protId].get("split")
                if protSet in self.remapSets:
                    self.proteins[protId]["split"] = self.remapSets[protSet]
            print "Remapped splits", Counter([x.get("split") for x in self.proteins.values()])
        else:
            loading.loadSplit(os.path.join(options.dataPath, "data"), proteins, self.allowMissing)
        if fold != None:
            makeFolds.loadFolds(proteins, os.path.join(options.dataPath, "folds", "training_folds_170125.tsv.gz" if task != "cafapi" else "CAFA_PI_training_folds_180417.tsv.gz"))
        loading.defineSets(proteins, cafaTargets, fold=fold, limitTrainingToAnnotated = task != "cafapi")
    
    def getTopTerms(self, counts, num=1000):
        return sorted(counts.items(), key=operator.itemgetter(1), reverse=True)[0:num]
        

class CAFAPITask(Task):
    def __init__(self):
        pass