import os
import gzip

def onError(message, errors):
    if errors == "strict":
        raise Exception(message)
    else:
        print "WARNING: " + message

def getFoldDirs(inPath, foldPattern, numFolds):
    return [(i, os.path.join(inPath, foldPattern.replace("{NUMBER}", str(i)))) for i in range(numFolds)]        

def collect(inPath, numFolds, foldPattern, errors):
    #predictions = {"devel":[], "test":[]}
    #seenIds = {"devel":set(), "test":set()}
    outFiles = {}
    for setName in ("devel", "test"):
        outFiles[setName] = gzip.open(os.path.join(inPath, setName + "allfolds-predicted.tsv.gz"))
    logLines = {i:{"cv":None, "command":None, "best_args":None, "avg_devel":None, "avg_test":None} for i in range(numFolds)}
    patterns = {"cv":"Cross-validation sets", "command":"Command line:", "best_args":""}
    for i, foldDir in getFoldDirs(inPath, foldPattern, numFolds):
        logPath = os.path.join(foldDir, "log.txt")
        if not os.path.exists(logPath):
            onError("Log file " + logPath + " not found", errors)
            continue
        with open(logPath, "rt") as f:
            for line in f:
                if line.
        
    for i, foldDir in getFoldDirs(inPath, foldPattern, numFolds):
        foldDir = os.path.join(inPath, foldPattern.replace("{NUMBER}", str(i)))
        print "Processing fold", i, foldDir
        if not os.path.exists(foldDir):
            onError("Directory " + foldDir + " not found", errors)
            continue
        for setName in ("devel", "test"):
            predPath = os.path.join(foldDir, setName + "-predictions.tsv.gz")
            print "Reading predictions for fold", i, "from", predPath
            if os.path.exists(predPath):
                with gzip.open(predPath, "rt") as f:
                    for line in f:
                        outFiles[setName].write(line)
            else:
                onError("Result file '" + predPath + "' not found")
        

def run(inPath, numFolds):
    validate(inPath, numFolds=numFolds)

if __name__=="__main__":       
    from optparse import OptionParser
    optparser = OptionParser(description="")
    optparser.add_option("-a", "--actions", default=None, help="One or more of 'build', 'classify' or 'statistics'")
    optparser.add_option("-i", "--input", default=None, help="")
    optparser.add_option("-p", "--dataPath", default=os.path.expanduser("~/data/CAFA3/data"), help="The main directory for the data files")
    optparser.add_option("-n", "--numFolds", default=10, type=int)
    optparser.add_option("-f", "--foldPattern", default="fold{NUMBER}")
    
    
    optparser.add_option("-f", "--features", default="all", help="Comma-separated list of feature group names. Use 'all' for all feature groups and '-name' to remove groups.")
    optparser.add_option("-l", "--limit", default=None, type=int, help="Limit the number of proteins to read.")
    optparser.add_option("-t", "--terms", default=100, type=int, help="The number of top most common GO terms to use as labels")
    optparser.add_option("-o", "--output", default=None, help="The output directory")
    optparser.add_option('-c','--classifier', help='', default="ensemble.RandomForestClassifier")
    optparser.add_option('-r','--args', help='', default="{'random_state':[1], 'n_estimators':[10], 'n_jobs':[1], 'verbose':[3]}")
    #optparser.add_option("--multioutputclassifier", default=False, action="store_true", help="Use the MultiOutputClassifier to train a separate classifier for each label")
    optparser.add_option("--singleLabelJobs", default=None, type=int, help="Number of jobs for SingleLabelClassification")
    optparser.add_option("--testSet", default=False, action="store_true", help="Classify the test set")
    optparser.add_option("--clear", default=False, action="store_true", help="Remove the output directory if it already exists")
    optparser.add_option("--targets", default="skip", help="How to include the CAFA target proteins, one of 'skip', 'overlap' or 'separate'")
    optparser.add_option("--negatives", default=False, action="store_true", help="Write negative predictions in the result files")
    (options, args) = optparser.parse_args()