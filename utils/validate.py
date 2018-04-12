import sys, os
from collections import defaultdict
import gzip
import csv

def loadCAFA(cafaDir, modelNumber, baselineConf, termTag, counts, tag):
    modelTag = "_" + str(modelNumber) + "_"
    proteins = {}
    for filename in os.listdir(cafaDir):
        filePath = os.path.join(cafaDir, filename)
        if modelTag in filePath:
            print "Reading", filePath
            counts["model-files-" + tag] += 1
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
                        counts["skipped-term-tag-" + tag + rowTermTag] += 1
                        continue
                    if cafaId not in proteins:
                        proteins[cafaId] = {}
                    if termId in proteins[cafaId]:
                        counts["duplicate-terms-" + tag] += 1
                    proteins[cafaId][termId] = confidence
                    counts["terms-total-" + tag] += 1
                    if confidence == baselineConf:
                        counts["terms-baseline-conf-" + tag] += 1
                    else:
                        counts["terms-non-baseline-conf-" + tag] += 1
    return proteins

def loadTSV(tsvPath, baselineConf, counts, tag):
    proteins = {}
    rowCount = 0
    with gzip.open(tsvPath, "rt") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            #print row, (baselineConf, baselineConf == row["confidence"])
            rowCount += 1
            if rowCount % 1000000 == 0:
                print "Processing row", rowCount
            if row["cafa_ids"] in ("", None):
                counts["skipped-terms-non-cafa-target-" + tag] += 1
                continue
            cafaIds = row["cafa_ids"].split(",")
            cafaIds = [x for x in cafaIds if not x.startswith("M")]
            if len(cafaIds) == 0:
                counts["moonlight-only-predictions-" + tag] += 1
                continue
            assert len(cafaIds) == 1, cafaIds
            cafaId = cafaIds[0]
            if row["predicted"] == "1":
                if cafaId not in proteins:
                    proteins[cafaId] = {}
                termId = row["label"]
                confidence = '%.2f' % float(row["confidence"])
                if termId in proteins[cafaId]:
                    counts["duplicate-terms-" + tag] += 1
                proteins[cafaId][termId] = confidence
                counts["terms-total-" + tag] += 1
                if confidence == baselineConf:
                    counts["terms-baseline-conf-" + tag] += 1
                else:
                    counts["terms-non-baseline-conf-" + tag] += 1
    return proteins

def loadSource(source, modelNumber, baselineConf, termTag, counts, tag):
    if os.path.isdir(source):
        return loadCAFA(source, modelNumber, baselineConf, termTag, counts, tag)
    else:
        return loadTSV(source, baselineConf, counts, tag)

def validate(sourceA, sourceB, modelNumber, baselineConf, termTag):
    counts = defaultdict(int)
    print "Reading A proteins from", sourceA
    proteinsA = loadSource(sourceA, modelNumber, baselineConf, termTag, counts, "A")
    print "Reading B proteins from", sourceB
    proteinsB = loadSource(sourceB, modelNumber, baselineConf, termTag, counts, "B")
    counts["proteins-A"] = len(proteinsA)
    counts["proteins-B"] = len(proteinsB)
    allProteinKeys = sorted(set(proteinsA.keys() + proteinsB.keys()))
    numProteinKeys = len(allProteinKeys)
    protCount = 0
    for protId in allProteinKeys:
        protCount += 1
        if protCount % 10000 == 0:
            print "Processing protein", protId, "(" + str(protCount) + "/" + str(numProteinKeys) + ")"
        if protId not in proteinsA:
            counts["protein-missing-A"] += 1
        if protId not in proteinsB:
            counts["protein-missing-B"] += 1
        protA = proteinsA.get(protId, {})
        protB = proteinsB.get(protId, {})
        for label in sorted(set(protA.keys() + protB.keys())):
            confA = protA.get(label)
            confB = protB.get(label)
            if confA != None and confB != None:
                counts["matching-terms-total"] += 1
                if confA != confB:
                    counts["matching-terms-different-conf"] += 1
                    #print key, label, (conf1, conf2)
                elif confA == baselineConf:
                    assert confB == baselineConf
                    counts["matching-terms-baseline-conf"] += 1
                else:
                    assert confA == confB
                    counts["matching-terms-non-baseline-conf"] += 1
            elif confA == None:
                counts["missing-terms-total-A"] += 1
                if confB == baselineConf:
                    counts["missing-baseline-conf-terms-A"] += 1
                else:
                    counts["missing-non-baseline-conf-terms-A"] += 1
            elif confB == None:
                counts["missing-terms-total-B"] += 1
                if confA == baselineConf:
                    counts["missing-baseline-conf-terms-B"] += 1
                else:
                    counts["missing-non-baseline-conf-terms-B"] += 1
    print "Statistics for model number", modelNumber
    print "Source A:", sourceA
    print "Source B:", sourceB
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