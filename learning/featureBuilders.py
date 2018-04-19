import gzip
import csv
import re
import os

###############################################################################
# Base Classes
###############################################################################

class FeatureBuilder:
    def __init__(self, debug=False):
        self.debug = debug
    
    def build(self, proteins):
        pass
    
    def setFeature(self, protein, feature, value=1):
        protein["features"][feature] = value
    
    def getMatchingPaths(self, inPaths, patterns):
        matching = []
        for inPath in inPaths:
            if not os.path.exists(inPath):
                raise Exception("Cannot find path " + str(inPath))
            for filename in sorted(os.listdir(inPath)):
                if filename.endswith("~"):
                    continue
                for pattern in patterns:
                    if pattern.match(filename) != None:
                        matching.append(os.path.join(inPath, filename))
                        break
        return matching
    
    def beginCoverage(self, proteins):
        self.numProteins = {"total":len(proteins)}
        for protein in proteins:
            protSets = ",".join(protein["sets"])
            if protSets not in self.numProteins:
                self.numProteins[protSets] = 0
            self.numProteins[protSets] += 1
        self.coverage = {"total":set()}
    
    def addToCoverage(self, protein):
        protId = protein["id"]
        self.coverage["total"].add(protId)
        protSets = ",".join(protein["sets"])
        if protSets not in self.coverage:
            self.coverage[protSets] = set()
        self.coverage[protSets].add(protId)
        #self.coveredIds.add(protId)
    
    def finishCoverage(self):
        #print self.__class__.__name__, "coverage =", float(len(self.coveredIds)) / self.numProteins, "included/total =" [len(self.coveredIds), self.numProteins]
        counts = {x:len(self.coverage[x]) for x in self.coverage}
        print self.__class__.__name__, "coverage (included, all, fraction) =", {x:(counts[x], self.numProteins[x], float(counts[x]) / self.numProteins[x]) for x in sorted(self.numProteins.keys())}
        self.numProteins = None
        self.coverage = None 

class MultiFileFeatureBuilder(FeatureBuilder):
    def __init__(self, inPaths, filePatterns, tag, message, debug=False):
        FeatureBuilder.__init__(self, debug=debug)
        self.tag = tag
        self.inPaths = inPaths
        self.filePatterns = filePatterns
        self.message = message
    
    def build(self, proteins):
        print self.message
        protById = {}
        for protein in proteins:
            protById[protein["id"]] = protein
        assert len(protById) == len(proteins)
        self.beginCoverage(proteins)
        for filePath in self.getMatchingPaths(self.inPaths, self.filePatterns):
            print "Reading", filePath
            try:
                self.buildForFile(filePath, protById)
            except IOError as e:
                print e
                print "Error reading file", filePath
                if not self.debug:
                    raise e
        self.finishCoverage()
    
    def buildForFile(self, filePath, protById):
        raise NotImplementedError()

class KeyValueFeatureBuilder(MultiFileFeatureBuilder):
    def __init__(self, inPaths, filePatterns, tag, message, skipHeader=True):
        MultiFileFeatureBuilder.__init__(self, inPaths, filePatterns, tag, message)
        self.skipHeader = skipHeader
    
    def buildForFile(self, filePath, protById):
        with gzip.open(filePath, "rt") as f:
            if self.skipHeader:
                f.readline() # Skip the headers
            for line in f:
                #print line.strip().split("\t")
                key, value = line.strip().split("\t")
                protein = protById.get(key)
                if protein is not None:
                    self.addToCoverage(protein)
                    self.setValue(protein, value)
    
    def setValue(self, protein, value):
        protein["features"][self.tag + ":value"] = float(value)

class CSVFeatureBuilder(MultiFileFeatureBuilder):
    def __init__(self, inPaths, filePatterns, tag, message, protColumn, columns=None, delimiter='\t', debug=False):
        MultiFileFeatureBuilder.__init__(self, inPaths, filePatterns, tag, message, debug=debug)
        self.columns = columns
        self.protColumn = protColumn
        self.delimiter = delimiter
    
    def buildForFile(self, filePath, protById):
        with gzip.open(filePath, "rt") as f:
            reader = csv.DictReader(f, delimiter=self.delimiter, fieldnames=self.columns)
            currentId = None
            found = False
            features = None
            for row in reader:
                if row[self.protColumn] != currentId:
                    currentId = row[self.protColumn]
                    found = currentId in protById
                    features = protById[currentId]["features"] if found else None
                    if found:
                        self.addToCoverage(protById[currentId])
                if found:
                    try:
                        self.setRow(features, row, filePath)
                    except ValueError as e:
                        print e
                        print "Error reading row from file", filePath
                        if not self.debug:
                            raise e
    
    def setRow(self, features, row, filePath):
        for key in row:
            if key != self.protColumn:
                features[self.tag + ":" + key] = float(row[key])
    
###############################################################################
# Feature Builders
###############################################################################

class BlastFeatureBuilder(CSVFeatureBuilder):
    def __init__(self, inPaths, tag="BLAST", debug=False): 
        filePatterns = (re.compile("target.[0-9]+.features_tsv\.gz"), re.compile("Swissprot\_sequence\_[0-9]\.features\_tsv\.gz"), re.compile("sequence\_[0-9]+\.fasta\.features\_tsv\.gz"))
        columns = ["Uniprot_ID query","Unknown_A","Unknown_B","Unknown_C","Matched Uniprot_ID","Matched Uniprot_ACC","Hsp_hit-len","Hsp_align-len","Hsp_bit-score","Hsp_score","Hsp_evalue","hsp.query_start","hsp.query_end","Hsp_hit-from","Hsp_hit-to","Hsp_query-frame","Hsp_hit-frame","Hsp_identity","Hsp_positives","Hsp_gaps"]
        CSVFeatureBuilder.__init__(self, inPaths, filePatterns, tag, "Building BLAST features", "Uniprot_ID query", columns, debug=debug)
    
    def setRow(self, features, row, filePath):
        features[self.tag + ":Hsp_score:" + row["Matched Uniprot_ID"]] = float(row["Hsp_score"])

class TaxonomyFeatureBuilder(KeyValueFeatureBuilder):
    def __init__(self, inPaths):
        filePatterns = [re.compile("map\_.+\_taxonomy\.tsv\.gz")]
        KeyValueFeatureBuilder.__init__(self, inPaths, filePatterns, "TAX", "Building taxonomy features", skipHeader=True)
    
    def setValue(self, protein, value):
        features = protein["features"]
        for taxonomyLevel in value.split(","):
            features[self.tag + ":" + taxonomyLevel] = 1

class NucPredFeatureBuilder(CSVFeatureBuilder):
    def __init__(self, inPaths):
        filePatterns = [re.compile("CAFA3\_nucPred\.tsv\.gz"), re.compile("training\_nucPred\.tsv\.gz")]
        CSVFeatureBuilder.__init__(self, inPaths, filePatterns, "NUC", "Building nucPred features", "Sequence-ID")

class PredGPIFeatureBuilder(CSVFeatureBuilder):
    def __init__(self, inPaths, debug=False):
        filePatterns = [re.compile("new\_CAFA3\_predGPI\.tsv\.gz"), re.compile("new\_training\_predGPI\.tsv\.gz")]
        CSVFeatureBuilder.__init__(self, inPaths, filePatterns, "GPI", "Building predGPI features", "protein_id", debug=debug)

class NetAcetFeatureBuilder(CSVFeatureBuilder):
    def __init__(self, inPaths):
        filePatterns = [re.compile(".+\_sequence\_NetAcet\_extract\.tsv\.gz")]
        CSVFeatureBuilder.__init__(self, inPaths, filePatterns, "NET", "Building NetAcet features", "id", ["id"] + ["col_" + str(x) for x in range(2,14)])

class NGramFeatureBuilder(CSVFeatureBuilder):
    def __init__(self, inPaths):
        filePatterns = [re.compile(".+-prot2ngram\.tsv\.gz")]
        CSVFeatureBuilder.__init__(self, inPaths, filePatterns, "NGRAM", "Building NGram features", "PROT_ID")
    
    def setRow(self, features, row, filePath):
        features[self.tag + ":NGRAM_ID=" + row["NGRAM_ID"]] = 1
    
class InterProScanFeatureBuilder(CSVFeatureBuilder):
    def __init__(self, inPaths, debug=False):
        filePatterns = [re.compile(".+_noGO.tsv.gz"), re.compile(".+_GO.tsv.gz")]
        CSVFeatureBuilder.__init__(self, inPaths, filePatterns, "IPS", "Building InterProScan features", "protein_id", debug=debug)
    
    def setRow(self, features, row, filePath):
        # Choose the primary identifier for the feature
        for idType in ("GOid", "ac", "db"):
            if idType in row:
                break
        # Choose the numeric value to use for the feature
        for valueType in ("score", "evalue", "bin"):
            if valueType in row:
                break
        # Use the chosen numeric value for the feature or use 1.0 if this is a binary feature
        value = 1.0 if valueType == "bin" else float(row[valueType])
        # Define the full name of the feature
        fileTag = os.path.basename(filePath).split(".")[0]
        name = self.tag + ":" + fileTag + ":" + idType + "=" + row[idType].replace(":", "_") + ":" + valueType
        # If the value is for a motif match, add this information to the feature name
        if "motifNumber" in row:
            name += ":motif=" + row["motifNumber"]
        # Add the feature for the protein
        features[name] = value

class FunTaxISFeatureBuilder(CSVFeatureBuilder):
    def __init__(self, inPaths):
        filePatterns = [re.compile(".+\_FunTaxIS\.tsv\.gz")]
        CSVFeatureBuilder.__init__(self, inPaths, filePatterns, "FUN", "Building FunTaxIS features", None)
        mapFilePatterns = [re.compile("map\_.+\_organism\.tsv\.gz")]
        self.mapping = self.readMapping(inPaths, mapFilePatterns)
    
    def buildForFile(self, filePath, protById):
        with gzip.open(filePath, "rt") as f:
            reader = csv.DictReader(f, delimiter=self.delimiter, fieldnames=self.columns)
            current_ncbitax_id = None
            currentProteins = None
            for row in reader:
                ncbitax_id = row["ncbitax_id"]
                if ncbitax_id != current_ncbitax_id:
                    current_ncbitax_id = ncbitax_id
                    currentProteins = []
                    if ncbitax_id in self.mapping:
                        for protId in self.mapping[ncbitax_id]:
                            if protId in protById:
                                self.addToCoverage(protById[protId])
                                currentProteins.append(protById[protId])
                    if len(currentProteins) == 0:
                        currentProteins = None
                if currentProteins != None:
                    baseName = self.tag + ":" + row["go_id"]
                    if row["conf"] == "None":
                        f1Name, f1Value = baseName + ":conf=None", 1
                    else:
                        f1Name, f1Value = baseName + ":conf", float(row["conf"])
                    f2Name, f2Value = baseName + ":no_protein", int(row["no_protein"])
                    for protein in currentProteins:
                        features = protein["features"]
                        features[f1Name] = f1Value
                        features[f2Name] = f2Value
    
    def readMapping(self, mapPaths, mapFilePatterns):
        print "Reading FunTaxIS Mapping"
        mapping = {}
        for filePath in self.getMatchingPaths(mapPaths, mapFilePatterns):
            print "Reading", filePath
            with gzip.open(filePath, "rt") as f:
                reader = csv.DictReader(f, delimiter=self.delimiter)
                for row in reader:
                    symbol, taxId = row["symbol"], row["ncbitax_id"]
                    #assert "_" in symbol, row
                    if taxId not in mapping:
                        mapping[taxId] = []
                    mapping[taxId].append(symbol)
        return mapping
    
class UniprotFeatureBuilder(FeatureBuilder):
    def __init__(self, inPath):
        FeatureBuilder.__init__(self)
        self.loadSimilar(inPath)
    
    def build(self, proteins):
        print "Building Uniprot similar.txt features"
        self.beginCoverage([x.get("id") for x in proteins])
        for protein in proteins:
            protId = protein["id"]
            if protId in self.data:
                self.addToCoverage(protein)
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