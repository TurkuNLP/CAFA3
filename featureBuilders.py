import gzip
import csv
import re
import os

class FeatureBuilder:
    def __init__(self):
        pass
    
    def build(self, proteins):
        pass
    
    def setFeature(self, protein, feature, value=1):
        protein["features"][feature] = value
    
    def getMatchingPaths(self, inPaths, patterns):
        matching = []
        for inPath in inPaths:
            for filename in sorted(os.listdir(inPath)):
                for pattern in patterns:
                    if pattern.match(filename) != None:
                        matching.append(os.path.join(inPath, filename))
                        break
        return matching
    
    def beginCoverage(self, protIds):
        self.numProteins = len(set(protIds))
        self.coveredIds = set()
    
    def addToCoverage(self, protId):
        self.coveredIds.add(protId)
    
    def finishCoverage(self):
        print self.__class__.__name__, "coverage =", float(len(self.coveredIds)) / self.numProteins, [len(self.coveredIds), self.numProteins]
        self.numProteins = None
        self.coveredIds = None 

class UniprotFeatureBuilder(FeatureBuilder):
    def __init__(self, inPath):
        self.loadSimilar(inPath)
    
    def build(self, proteins):
        print "Building Uniprot similar.txt features"
        self.beginCoverage([x.get("id") for x in proteins])
        for protein in proteins:
            protId = protein["id"]
            if protId in self.data:
                self.addToCoverage(protId)
                for section in ("sub", "fam"):
                    for feature in self.data[protId][section]:
                        self.setFeature(protein, "SIM:" + section + ":" + feature, 1)
        self.finishCoverage()
    
    def loadSimilar(self, inPath):
        self.data = {}
        print "Loading uniprot similar.txt"
        with open(inPath, "rt") as f:
            section = None
            group = None
            subgroup = ""
            #lineCount = 0
            for line in f:
                if line.startswith("I. Domains, repeats and zinc fingers"):
                    section = "sub"
                elif line.startswith("II. Families"):
                    section = "fam"
                elif section == None: # Not yet in the actual data
                    continue
                elif line.startswith("----------"): # end of file
                    break
                elif line.strip() == "":
                    continue
                #elif line.strip() == "": # empty line ends block
                #    group = None
                elif line[0].strip() != "": # a new group (family or domain)
                    group = line.strip().replace(" ", "-").replace(",",";")
                    subgroup = ""
                elif line.startswith("  ") and line[2] != " ":
                    subgroup = "<" + line.strip().replace(" ", "-").replace(",",";") + ">"
                elif line.startswith("     "):
                    assert group != None, line
                    items = [x.strip() for x in line.strip().split(",")]
                    items = [x for x in items if x != ""]
                    protIds = [x.split()[0] for x in items]
                    for protId in protIds:
                        if protId not in self.data:
                            self.data[protId] = {"sub":set(), "fam":set()}
                        self.data[protId][section].add(group + subgroup)

class BlastFeatureBuilder(FeatureBuilder):
    def __init__(self, inPaths, tag="BLAST"):
        self.tag = tag
        self.inPaths = inPaths
        self.filePatterns = (re.compile("target.[0-9]+.features_tsv.gz"), re.compile("Swissprot\_sequence\_[0-9].features\_tsv.gz"))
        self.columns = ["Uniprot_ID query","Unknown_A","Unknown_B","Unknown_C","Matched Uniprot_ID","Matched Uniprot_ACC","Hsp_hit-len","Hsp_align-len","Hsp_bit-score","Hsp_score","Hsp_evalue","hsp.query_start","hsp.query_end","Hsp_hit-from","Hsp_hit-to","Hsp_query-frame","Hsp_hit-frame","Hsp_identity","Hsp_positives","Hsp_gaps"]
    
    def build(self, proteins):
        print "Building BLASTP features"
        protById = {}
        for protein in proteins:
            protById[protein["id"]] = protein
        self.beginCoverage(protById.keys())
        for filePath in self.getMatchingPaths(self.inPaths, self.filePatterns):
            print "Reading", filePath
            with gzip.open(filePath, "rt") as f:
                reader = csv.DictReader(f, delimiter='\t', fieldnames=self.columns)
                current = None
                found = False
                features = None
                for row in reader:
                    if row["Uniprot_ID query"] != current:
                        current = row["Uniprot_ID query"]
                        found = current in protById
                        features = protById[current]["features"] if found else None
                        if found:
                            self.addToCoverage(current)
                    if found:
                        features[self.tag + ":Hsp_score:" + row["Matched Uniprot_ID"]] = float(row["Hsp_score"])
        self.finishCoverage()

class MultiFileFeatureBuilder(FeatureBuilder):
    def __init__(self, inPaths, filePatterns, tag):
        self.tag = tag
        self.inPaths = inPaths
        self.filePatterns = filePatterns
    
    def build(self, proteins):
        print "Building InterProScan features"
        protById = {}
        for protein in proteins:
            protById[protein["id"]] = protein
        self.beginCoverage(protById.keys())
        for filePath in self.getMatchingPaths(self.inPaths, self.filePatterns):
            print "Reading", filePath
            with gzip.open(filePath, "rt") as f:
                reader = csv.DictReader(f, delimiter='\t')
                current = None
                found = False
                features = None
                for row in reader:
                    if row["Uniprot_ID query"] != current:
                        current = row["Uniprot_ID query"]
                        found = current in protById
                        features = protById[current]["features"] if found else None
                        if found:
                            self.addToCoverage(current)
                    if found:
                        features[self.tag + ":Hsp_score:" + row["Matched Uniprot_ID"]] = float(row["Hsp_score"])
        self.finishCoverage()

class KeyValueFeatureBuilder(FeatureBuilder):
    def __init__(self):
        self.skipHeader = None
        self.message = "Building features"
    
    def build(self, proteins):
        print self.message
        protById = {}
        for protein in proteins:
            protById[protein["id"]] = protein
        self.beginCoverage(protById.keys())
        for filePath in self.getMatchingPaths([self.inPath], [self.filePattern]):
            print "Reading", filePath
            with gzip.open(filePath, "rt") as f:
                if self.skipHeader:
                    f.readline() # Skip the headers
                for line in f:
                    #print line.strip().split("\t")
                    key, value = line.strip().split("\t")
                    protein = protById.get(key)
                    if protein is not None:
                        self.addToCoverage(protein["id"])
                        self.setValue(protein, value)
        self.finishCoverage()
    
    def setValue(self, protein, value):
        protein["features"][self.tag] = float(value)

class TaxonomyFeatureBuilder(KeyValueFeatureBuilder):
    def __init__(self, inPath):
        super(TaxonomyFeatureBuilder, self).__init__()
        self.message = "Building taxonomy features"
        self.inPath = inPath
        self.filePattern = re.compile("map\_.+\_taxonomy\.tsv\.gz") #"*_taxonomy_lineage.tsv.gz"
    
    def setValue(self, protein, value):
        features = protein["features"]
        for taxonomyLevel in value.split(","):
            features[self.tag + ":" + taxonomyLevel] = 1

class InterproScanFeatureBuilder(FeatureBuilder):
    def __init__(self, inPaths, tag="IPS"):
        self.tag = tag
        self.inPaths = inPaths
        self.filePatterns = [re.compile(".+_noGO.tsv.gz")]
    
    def build(self, proteins):
        print "Building InterProScan features"
        protById = {}
        for protein in proteins:
            protById[protein["id"]] = protein
        self.beginCoverage(protById.keys())
        for filePath in self.getMatchingPaths(self.inPaths, self.filePatterns):
            print "Reading", filePath
            with gzip.open(filePath, "rt") as f:
                reader = csv.DictReader(f, delimiter='\t')
                current = None
                found = False
                features = None
                for row in reader:
                    if row["Uniprot_ID query"] != current:
                        current = row["Uniprot_ID query"]
                        found = current in protById
                        features = protById[current]["features"] if found else None
                        if found:
                            self.addToCoverage(current)
                    if found:
                        features[self.tag + ":Hsp_score:" + row["Matched Uniprot_ID"]] = float(row["Hsp_score"])
        self.finishCoverage()