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
    
    def getMatchingPaths(self, inPath, pattern):
        matching = []
        for filename in sorted(os.listdir(self.inPath)):
            if self.filePattern.match(filename) != None:
                matching.append(os.path.join(self.inPath, filename))
        return matching

class UniprotFeatureBuilder(FeatureBuilder):
    def __init__(self, inPath):
        self.loadSimilar(inPath)
    
    def build(self, proteins):
        print "Building Uniprot similar.txt features"
        for protein in proteins:
            protId = protein["id"]
            if protId in self.data:
                for section in ("sub", "fam"):
                    for feature in self.data[protId][section]:
                        self.setFeature(protein, "SIM:" + section + ":" + feature, 1)
    
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
    def __init__(self, inPath):
        self.inPath = inPath
        self.filePattern = re.compile("target.[0-9]+.features_tsv.gz")
        self.columns = ["Uniprot_ID query","Unknown_A","Unknown_B","Unknown_C","Matched Uniprot_ID","Matched Uniprot_ACC","Hsp_hit-len","Hsp_align-len","Hsp_bit-score","Hsp_score","Hsp_evalue","hsp.query_start","hsp.query_end","Hsp_hit-from","Hsp_hit-to","Hsp_query-frame","Hsp_hit-frame","Hsp_identity","Hsp_positives","Hsp_gaps"]
    
    def build(self, proteins):
        print "Building BLASTP features"
        protById = {}
        for protein in proteins:
            protById[protein["id"]] = protein
        for filePath in self.getMatchingPaths(self.inPath, self.filePattern):
            print "Reading", os.path.basename(filePath)
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
                        features["BLAST:Hsp_score:" + row["Matched Uniprot_ID"]] = float(row["Hsp_score"])

class TaxonomyFeatureBuilder(FeatureBuilder):
    def __init__(self, inPath):
        self.inPath = inPath
        self.filePattern = re.compile("map\_.+\_taxonomy\.tsv\.gz") #"*_taxonomy_lineage.tsv.gz"
    
    def build(self, proteins):
        print "Building taxonomy features"
        protById = {}
        for protein in proteins:
            protById[protein["id"]] = protein
        for filePath in self.getMatchingPaths(self.inPath, self.filePattern):
            print "Reading", os.path.basename(filePath)
            with gzip.open(filePath, "rt") as f:
                f.readline() # Skip the headers
                for line in f:
                    #print line.strip().split("\t")
                    symbol, taxonomy = line.strip().split("\t")
                    protein = protById.get(symbol)
                    if protein is not None:
                        features = protein["features"]
                        for level in taxonomy.split(","):
                            features["TAX:" + level] = 1
                        