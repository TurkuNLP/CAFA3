import sys, os
from collections import defaultdict
import gzip
import csv

def loadCAFA(cafaDir, modelNumber, baselineConf, counts):
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

def validate(cafaDir, tsvPath, modelNumber, baselineConf):
    counts = defaultdict(int)
    print "Reading CAFA proteins"
    proteins1 = loadCAFA(cafaDir, modelNumber, baselineConf, counts)
    print "Readgind TSV proteins"
    proteins2 = loadTSV(tsvPath, baselineConf, counts)
    allProteinKeys = sorted(set(proteins1.keys() + proteins2.keys()))
    numProteinKeys = allProteinKeys
    protCount = 0
    for key in allProteinKeys:
        protCount += 1
        print "Processing protein", key, "(" + str(protCount) + "/" + str(numProteinKeys) + ")"
        if key not in proteins1:
            counts["cafa-missing-protein"] += 1
        if key not in proteins2:
            counts["tsv-missing-protein"] += 1
        prot1 = proteins1.get(key, {})
        prot2 = proteins2.get(key, {})
        for label in sorted(set(prot1.keys() + prot2.keys())):
            conf1 = prot1.get(label)
            conf2 = prot2.get(label)
            if conf1 != None and conf1 == conf2:
                if conf1 == baselineConf:
                    counts["matching-baseline-terms"] += 1
                else:
                    counts["matching-non-baseline-terms"] += 1
            elif conf1 == None:
                counts["cafa-missing-label"] += 1
                if conf2 == baselineConf:
                    counts["cafa-missing-baseline-terms"] += 1
                else:
                    counts["cafa-missing-non-baseline-terms"] += 1
            elif conf2 == None:
                counts["tsv-missing-label"] += 1
                if conf1 == baselineConf:
                    counts["cafa-missing-baseline-terms"] += 1
                else:
                    counts["cafa-missing-non-baseline-terms"] += 1
    print "Statistics"
    for key in sorted(counts.keys()):
        print key + ": " + str(counts[key])

if __name__=="__main__":       
    from optparse import OptionParser
    optparser = OptionParser(description="")
    optparser.add_option("-a", "--cafaDir", default=None, help="")
    optparser.add_option("-b", "--tsvPath", default=None, help="")
    optparser.add_option("-m", "--model", default=None, type=int, help="")
    optparser.add_option("-c", "--baselineConf", default="1.01", help="")
    (options, args) = optparser.parse_args()
    
    validate(options.cafaDir, options.tsvPath, options.model, options.baselineConf)