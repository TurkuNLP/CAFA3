class FeatureBuilder:
    def __init__(self):
        pass
    
    def build(self, proteins):
        pass
    
    def setFeature(self, protein, feature, value=1):
        protein["features"][feature] = value

class UniprotFeatureBuilder(FeatureBuilder):
    def __init__(self, inPath):
        self.data = self.loadSimilar(inPath)
    
    def build(self, proteins):
        for protein in proteins:
            protId = protein["id"]
            if id in self.data:
                for section in ("sub", "fam"):
                    for feature in self.data[protId][section]:
                        self.setFeature(protein, feature, 1)
    
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