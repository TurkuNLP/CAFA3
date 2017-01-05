import csv

def convert(inPath, outPath):
    assert inPath != None and outPath != None
    assert inPath != outPath
    with open(inPath, "rt") as f:
        terms = []
        term = None
        for line in f:
            if line.startswith("["):
                term = None
                if line.startswith("[Term]"):
                    term = {}
                    terms.append(term)
            elif term != None and ":" in line:
                line = line.strip()
                tag, content = line.split(":", 1)
                term[tag] = content
                if tag == "namespace":
                    term["ns"] = "".join([x[0] for x in content.split("_")])
    with open(outPath, "wt") as f:
        dw = csv.DictWriter(f, ["id", "ns", "name"], delimiter='\t')
        dw.writeheader()
        dw.writerows(terms)   

if __name__=="__main__":       
    from optparse import OptionParser
    optparser = OptionParser(description="")
    optparser.add_option("-i", "--input", default=None, help="")
    optparser.add_option("-o", "--output", default=None, help="")
    (options, args) = optparser.parse_args()
    
    convert(options.input, options.output)