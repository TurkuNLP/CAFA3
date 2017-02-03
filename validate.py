import sys, os
from collections import defaultdict
import gzip
import csv

def loadCAFA(cafaDir, modelNumber, baselineConf, termTag, counts):
    modelTag = "_" + str(modelNumber) + "_"
    proteins = {}
    for filename in os.listdir(cafaDir):
        filePath = os.path.join(cafaDir, filename)
        if modelTag in filePath:
            print "Reading", filePath
            counts["model-files-cafa"] += 1
            with open(os.path.join(filePath), "rt") as f:
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
                    rowTermTag = termId.split(":")[0]
                    if rowTermTag != termTag:
                        counts["skipped-term-tag-cafa-" + rowTermTag] += 1
                        continue
                    if cafaId not in proteins:
                        proteins[cafaId] = {}
                    if termId in proteins[cafaId]:
                        counts["duplicate-terms-cafa"] += 1
                    proteins[cafaId][termId] = confidence
                    counts["terms-total-cafa"] += 1
                    if confidence == baselineConf:
                        counts["terms-baseline-conf-cafa"] += 1
                    else:
                        counts["terms-non-baseline-conf-cafa"] += 1
    return proteins

def loadTSV(tsvPath, baselineConf, counts):
    proteins = {}
    rowCount = 0
    with gzip.open(tsvPath, "rt") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            #print row, (baselineConf, baselineConf == row["confidence"])
            rowCount += 1
            if rowCount % 1000000 == 0:
                print "Processing row", rowCount
            if row["cafa_ids"] in ("", None):
                counts["proteins-non-cafa-target-tsv"] += 1
            cafaIds = row["cafa_ids"].split(",")
            cafaIds = [x for x in cafaIds if not x.startswith("M")]
            if len(cafaIds) == 0:
                counts["moonlight-only-predictions-tsv"] += 1
                continue
            assert len(cafaIds) == 1, cafaIds
            cafaId = cafaIds[0]
            if row["predicted"] == "1":
                if cafaId not in proteins:
                    proteins[cafaId] = {}
                termId = row["label"]
                confidence = '%.2f' % float(row["confidence"])
                if termId in proteins[cafaId]:
                    counts["duplicate-terms-tsv"] += 1
                proteins[cafaId][termId] = confidence
                counts["terms-total-tsv"] += 1
                if confidence == baselineConf:
                    counts["terms-baseline-conf-tsv"] += 1
                else:
                    counts["terms-non-baseline-conf-tsv"] += 1
    return proteins

def validate(cafaDir, tsvPath, modelNumber, baselineConf, termTag):
    counts = defaultdict(int)
    print "Reading CAFA proteins"
    proteins1 = loadCAFA(cafaDir, modelNumber, baselineConf, termTag, counts)
    print "Reading TSV proteins"
    proteins2 = loadTSV(tsvPath, baselineConf, counts)
    counts["proteins-cafa"] = len(proteins1)
    counts["proteins-tsv"] = len(proteins2)
    allProteinKeys = sorted(set(proteins1.keys() + proteins2.keys()))
    numProteinKeys = len(allProteinKeys)
    protCount = 0
    for key in allProteinKeys:
        protCount += 1
        if protCount % 10000 == 0:
            print "Processing protein", key, "(" + str(protCount) + "/" + str(numProteinKeys) + ")"
        if key not in proteins1:
            counts["protein-missing-cafa"] += 1
        if key not in proteins2:
            counts["protein-missing-tsv"] += 1
        prot1 = proteins1.get(key, {})
        prot2 = proteins2.get(key, {})
        for label in sorted(set(prot1.keys() + prot2.keys())):
            conf1 = prot1.get(label)
            conf2 = prot2.get(label)
            if conf1 != None and conf1 != None:
                counts["matching-terms-total"] += 1
                if conf1 != conf2:
                    counts["matching-terms-different-conf"] += 1
                    #print key, label, (conf1, conf2)
                elif conf1 == baselineConf:
                    counts["matching-terms-baseline-conf"] += 1
                else:
                    counts["matching-terms-non-baseline-conf"] += 1
            elif conf1 == None:
                counts["missing-terms-total-cafa"] += 1
                if conf2 == baselineConf:
                    counts["missing-baseline-conf-terms-cafa"] += 1
                else:
                    counts["missing-non-baseline-conf-terms-cafa"] += 1
            elif conf2 == None:
                counts["missing-terms-total-tsv"] += 1
                if conf1 == baselineConf:
                    counts["missing-baseline-conf-terms-tsv"] += 1
                else:
                    counts["missing-non-baseline-conf-terms-tsv"] += 1
    print "Statistics for model number", modelNumber
    print "CAFA model directory:", cafaDir
    print "TSV file:", tsvPath
    for key in sorted(counts.keys()):
        print key + ": " + str(counts[key])

if __name__=="__main__":       
    from optparse import OptionParser
    optparser = OptionParser(description="")
    optparser.add_option("-a", "--cafaDir", default=None, help="")
    optparser.add_option("-b", "--tsvPath", default=None, help="")
    optparser.add_option("-m", "--model", default=None, type=int, help="")
    optparser.add_option("-c", "--baselineConf", default="0.01", help="")
    optparser.add_option("-t", "--termTag", default="GO", help="")
    (options, args) = optparser.parse_args()
    
    validate(options.cafaDir, options.tsvPath, options.model, options.baselineConf, options.termTag)