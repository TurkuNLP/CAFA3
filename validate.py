import sys, os
from collections import defaultdict
import gzip
import csv

def loadCAFA(cafaDir, modelNumber, baselineConf):
    counts = defaultdict(int)
    counts["cafa_duplicate_terms"] = 0
    modelTag = "_" + str(modelNumber) + "_"
    proteins = {}
    for filename in os.listdir(cafaDir):
        if modelTag in filename:
            print "Reading", filename
            counts["cafa-model-files"] += 1
            with open(os.path.join(filename), "rt") as f:
                line = f.readline()
                assert line.startswith("AUTHOR")
                line = f.readline()
                assert line.startswith("MODEL " + str(modelNumber))
                line = f.readline()
                assert line.startswith("KEYWORDS")
                for line in f:
                    if line.startswith("END"):
                        continue
                    cafaId, termId, confidence = line.strip().split()
                    if cafaId not in proteins:
                        proteins[cafaId] = {}
                    if termId in proteins[cafaId]:
                        counts["cafa_duplicate_terms"] += 1
                    proteins[cafaId][termId] = confidence
                    counts["cafa-terms"] += 1
                    if confidence == baselineConf:
                        counts["cafa-baseline-terms"] += 1
                    else:
                        counts["cafa-non-baseline-terms"] += 1

def loadTSV(tsvPath, baselineConf, counts):
    proteins = {}
    with gzip.open(tsvPath, "rt") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            cafaIds = row["cafa_ids"].split(",")
            cafaIds = [x for x in cafaIds if not x.startswith("M")]
            assert len(cafaIds) == 1
            cafaId = cafaIds[0]
            if row["predicted"] == "1":
                if cafaId not in proteins:
                    proteins[cafaId] = {}
                termId = row["label"]
                confidence = row["confidence"]
                if termId in proteins[cafaId]:
                    counts["tsv_duplicate_terms"] += 1
                proteins[cafaId][termId] = confidence
                if confidence == baselineConf:
                    counts["tsv-baseline-terms"] += 1
                else:
                    counts["tsv-non-baseline-terms"] += 1
    return proteins
        

def validate(cafaDir, tsvPath):
    counts = d

if __name__=="__main__":       
    from optparse import OptionParser
    optparser = OptionParser(description="")
    optparser.add_option("-a", "--cafaDir", default=None, help="")
    optparser.add_option("-b", "--tsvPath", default=None, help="")
    optparser.add_option("-m", "--model", default=None, type=int, help="")
    optparser.add_option("-c", "--baselineConf", default="1.01", type=int, help="")
    (options, args) = optparser.parse_args()
    
    convert(options.input, options.output)