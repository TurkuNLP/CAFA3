import os
import gzip
import operator
import filecmp
import shutil

HEADER_LINE = "id\tlabel_index\tlabel\tpredicted\tconfidence\tgold\tmatch\tcafa_ids\r\n"

def onError(message, errors):
    if errors == "strict":
        raise Exception(message)
    else:
        print "WARNING: " + message

def getFoldDirs(inPath, foldPattern, numFolds):
    return [(i, os.path.join(inPath, foldPattern.replace("{NUMBER}", str(i)))) for i in range(numFolds)]        

def readLogs(inPath, foldPattern, numFolds, errors):
    print "Reading log files"
    patterns = ["Cross-validation sets", "Command line:", "Best classifier arguments:", "Best development set results:"]
    logLines = {i:{pattern:None for pattern in patterns + ["Test Average"]} for i in range(numFolds)} 
    for i, foldDir in getFoldDirs(inPath, foldPattern, numFolds):
        logPath = os.path.join(foldDir, "log.txt")
        if not os.path.exists(logPath):
            onError("Log file " + logPath + " not found", errors)
            continue
        with open(logPath, "rt") as f:
            for line in f:
                for pattern in patterns:
                    if pattern in line:
                        assert logLines[i][pattern] == None
                        logLines[i][pattern] = line.strip()
                if logLines[i]["Best development set results:"] != None and "Average:" in line and logLines[i]["Test Average"] == None:
                    logLines[i]["Test Average"] = line.strip()
    s = ""
    for pattern in patterns + ["Test Average"]:
        s += "*** " + pattern + " ***\n"
        for i in range(numFolds):
            s += str(i) + ":\t" + (logLines[i][pattern] if logLines[i][pattern] != None else "NOT FOUND") + "\n"
    argFolds = {}
    for i in range(numFolds):
        argLine = logLines[i]["Best classifier arguments:"]
        if argLine != None:
            argString = argLine.split("\t")[-1].strip()
            assert argString[0] == "B", argString
            if argString not in argFolds:
                argFolds[argString] = []
            argFolds[argString].append(i)
    argFolds = sorted([(len(argFolds[key]), argFolds[key], key) for key in argFolds], reverse=True)
    s += "Best arguments frequency\n"
    for argFold in argFolds:
        s += str(argFold) + "\n"
    print "Most common arguments", argFolds[0]
    return s, argFolds[0][1]

def checkHeaders(f, outFile):
    line = f.readline()
    assert line == HEADER_LINE, (HEADER_LINE,  line)
    if outFile != None:
        outFile.write(line)

def collect(inPath, numFolds, foldPattern, errors):
    #predictions = {"devel":[], "test":[]}
    #seenIds = {"devel":set(), "test":set()}
    logText, mostCommonArgsFolds = readLogs(inPath, foldPattern, numFolds, errors)
    print "Most common arguments are for folds", mostCommonArgsFolds
    with open(os.path.join(inPath, "logs.txt"), "wt") as f:
        f.write(logText)
    print "Merging predictions"
    outFiles = {}
    for setName in ("devel", "test"):
        outFiles[setName] = gzip.open(os.path.join(inPath, setName + "-allfolds-predicted.tsv.gz"), "wt")        
    chosenCAFAPath = None
    foldDirs = getFoldDirs(inPath, foldPattern, numFolds)
    for i, foldDir in foldDirs:
        foldDir = os.path.join(inPath, foldPattern.replace("{NUMBER}", str(i)))
        print "Processing fold", i, foldDir
        if not os.path.exists(foldDir):
            onError("Directory " + foldDir + " not found", errors)
            continue
        for setName in ("devel", "test"):
            predPath = os.path.join(foldDir, setName + "-predictions.tsv.gz")
            print "Reading", setName, "predictions for fold", i, "from", predPath
            if os.path.exists(predPath):
                with gzip.open(predPath, "rt") as f:
                    checkHeaders(f, outFiles[setName] if (i == foldDirs[0][0]) else None)
                    for line in f:
                        outFiles[setName].write(line)
            else:
                onError("Result file '" + predPath + "' not found")
        if i in mostCommonArgsFolds:
            foldCAFAPath = os.path.join(foldDir, "cafa-predictions.tsv.gz")
            if os.path.exists(foldCAFAPath):
                if chosenCAFAPath != None and os.path.exists(chosenCAFAPath):
                    sizes = [os.path.getsize(x) for x in (foldCAFAPath, chosenCAFAPath)]
                    print "Comparing", (foldCAFAPath, chosenCAFAPath), sizes
                    assert sizes[0] == sizes[1]
                else:
                    chosenCAFAPath = foldCAFAPath
                    print "Reading CAFA predictions for fold", i, "from", chosenCAFAPath
                    with gzip.open(chosenCAFAPath, "rt") as f:
                        checkHeaders(f, None) # skip headers
                        for line in f:
                            outFiles["test"].write(line)
            else:
                onError("Result file '" + foldCAFAPath + "' not found")
    for outFile in outFiles.values():
        outFile.close()

if __name__=="__main__":       
    from optparse import OptionParser
    optparser = OptionParser(description="")
    optparser.add_option("-i", "--input", default=None, help="")
    optparser.add_option("-p", "--dataPath", default=os.path.expanduser("~/data/CAFA3/data"), help="The main directory for the data files")
    optparser.add_option("-n", "--numFolds", default=10, type=int)
    optparser.add_option("-f", "--foldPattern", default="fold{NUMBER}")
    optparser.add_option("-e", "--errors", default="strict")
    (options, args) = optparser.parse_args()
    
    collect(options.input, options.numFolds, options.foldPattern, options.errors)