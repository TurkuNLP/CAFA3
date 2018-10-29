import os
from collections import defaultdict
import learning.evaluateFile as evaluateFile
import learning.loading as loading
import learning.evaluation as evaluation
import shutil
from utils import Stream
import itertools
from task.tasks import Task

def clearKeys(proteins, keys):
    for protId in proteins:
        protein = proteins[protId]
        for key in keys:
            if key in protein:
                del protein[key]

def getCombinations(items):
    combinations = []
    for i in xrange(1, len(items) + 1):
        els = [list(x) for x in itertools.combinations(items, i)]
        combinations.extend(els)
    return combinations

def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)

def combineConf(protein, labels, predKeys, combKey):
    combConfKey = combKey + "_conf"
    assert combConfKey not in protein
    protein[combConfKey] = {}
    combConfs = protein[combConfKey]
    assert combKey + "_sources" not in protein
    protein[combKey + "_sources"] = {}
    combSources = protein[combKey + "_sources"]
    for key in predKeys:
        preds = protein.get(key, {})
        predConfs = protein.get(key + "_conf", {})
        for label in labels:
            predConf = predConfs.get(label)
            if predConf != None:
                if label not in combConfs:
                    combConfs[label] = []
                combConfs[label].append(predConf)
            if label in preds:
                if label not in combSources:
                    combSources[label] = []
                combSources[label].append(key)
    protein[combConfKey] = {x:mean(combConfs[x]) for x in combConfs}
                        
def combinePred(proteins, predKeys, combKey, mode="AND", limitToSets=None):
    assert mode in ("AND", "OR", "SINGLE")
    if mode == "SINGLE":
        assert len(predKeys) == 1
    counts = defaultdict(int)
    for protId in proteins:
        protein = proteins[protId]
        if limitToSets != None and not any(x in limitToSets for x in protein["sets"]):
            counts["out-of-sets"] += 1
            continue
        counts["proteins"] += 1
        assert combKey not in protein
        protein[combKey] = {}
        if mode == "OR" or mode == "SINGLE":
            for key in predKeys:
                if key in protein:
                    counts["predictions-mode-" + mode] += 1
                    protein[combKey].update(protein[key])
                else:
                    counts["no-prediction-for-" + key] += 1
        elif mode == "AND":
            missing = False
            predLabelSets = []
            for key in predKeys:
                if key not in protein:
                    missing = True
                    counts["no-prediction-for-" + key] += 1
                else:
                    predLabelSets.append(set(protein[key].keys()))
            if not missing:
                counts["predictions-mode-" + mode] += 1
                protein[combKey] = {x:1 for x in set.intersection(*predLabelSets)} #{x:1 for x in pred1.intersection(pred2)}
        combineConf(protein, sorted(protein[combKey].keys()), predKeys, combKey)            
    print "Combined predictions, mode =", mode, "counts =", dict(counts)

def loadPredictionFiles(proteins, predPath, predKey, useCafa, task, predTags, readGold):
    inPaths = []
    if os.path.isfile(predPath):
        inPaths = [predPath]
    else:
        for setName in (("devel", "test", "cafa") if useCafa else ("devel", "test")):
            if setName == "test" and task == "cafapi":
                continue
            candidateNames = [setName + predTag + ".tsv.gz" for predTag in predTags[setName]]
            found = None
            for candidateName in candidateNames:
                inPath = os.path.join(predPath, candidateName)
                if os.path.exists(inPath):
                    found = inPath
                    break
            if found == None:
                raise Exception("No input file for set '" + setName + "' at " + predPath + ". Input file must be one of " + candidateNames)
            inPaths.append(found)
    for inPath in inPaths:
        evaluateFile.loadPredictions(proteins, inPath, limitToSets=["devel","test","cafa"] if useCafa else ["devel","test"], readGold=readGold, predKey=predKey, confKey=predKey + "_conf")
    
def combine(dataPath, inputs, outDir=None, cafaTargets="external", modes="AND,OR", baselineCutoff=1, numTerms=5000, clear=False, useOutFiles=True, taskName="cafa3", debug=False):
    if isinstance(inputs, basestring):
        inputs = [x.strip() for x in inputs.split(",")]
    for i in range(len(inputs)):
        if isinstance(inputs[i], basestring):
            item = [x.strip() for x in inputs[i].split(":")]
            assert len(item) == 3
            item = {"name":item[0], "type":item[1], "path":item[2]}
            inputs[i] = item
    if isinstance(modes, basestring):
        modes = [x.strip() for x in modes.split(",")]
    for mode in modes:
        assert mode in ("AND", "OR")
    
    if outDir != None:
        if clear and os.path.exists(outDir):
            print "Removing output directory", outDir
            shutil.rmtree(outDir)
        if not os.path.exists(outDir):
            print "Making output directory", outDir
            os.makedirs(outDir)
        Stream.openLog(os.path.join(outDir, "log.txt"))
    
    print "==========", "Ensemble", "=========="
    
    # Define the task
    task = Task.getTask(taskName)
    task.setDataPath(dataPath)
    task.setDebug(debug)
    if numTerms != None:
        task.numTerms = numTerms
    print "Task:", taskName
    #exampleFilePath = os.path.join(outDir, "examples.json.gz")
    task.loadProteins(cafaTargets)
    task.loadSplit()
    limitTerms = set([x[0] for x in task.topTerms])
    
    predKeys = [x["name"] for x in inputs]
    for item in inputs:
        print "Loading input", item
        assert item["type"] in ("cls", "bl")
        if (item["type"] == "cls"):
            loadPredictionFiles(task.proteins, item["path"], item["name"], cafaTargets != "skip", task, {"cafa":("_targets", "-predictions"), "devel":("_pred", "-predictions"), "test":("_pred", "-predictions")}, readGold=item == inputs[0])
        else:
            assert baselineCutoff >= 1
            loading.loadBaseline(item["path"], task.proteins, item["name"], baselineCutoff, limitTerms=limitTerms, useCafa=cafaTargets != "skip")
            
#     if nnInput != None:
#         print "Loading neural network predictions from", nnInput
#         loadPredictionFiles(task.proteins, nnInput, "nn_pred", cafaTargets != "skip", task, {"cafa":"_targets", "devel":"_pred", "test":"_pred"}, readGold=False)
#         predKeys += ["nn_pred"]
#     if clsInput != None:
#         print "Loading classifier predictions"
#         loadPredictionFiles(task.proteins, clsInput, "cls_pred", cafaTargets != "skip", task, {"cafa":"-predictions", "devel":"-predictions", "test":"-predictions"}, readGold=True)
#         predKeys += ["cls_pred"]
#     if baselineInput != None:
#         print "Loading baseline predictions from", baselineInput
#         assert baselineCutoff >= 1
#         loading.loadBaseline(baselineInput, task.proteins, "baseline_pred", baselineCutoff, limitTerms=limitTerms, useCafa=cafaTargets != "skip")
#         predKeys += ["baseline_pred"]
    
    coverage = {x:0 for x in predKeys}
    coverage["total"] = len(task.proteins)
    for protId in task.proteins:
        for predKey in predKeys:
            if predKey in task.proteins[protId]:
                coverage[predKey] += 1
    print "Coverage:", coverage
    

    print "===============", "Combining predictions", "===============" 
    combKey = "comb_pred"
    combConfKey = "comb_pred_conf"
    combinations = getCombinations(predKeys)
    numCombinations = len(combinations)
    print "Testing", numCombinations, "combinations"
    for i in range(len(combinations)):
        print
        print "******************", "Combination", str(i + 1) + "/" + str(numCombinations), combinations[i], "******************"
        for mode in (modes if len(combinations[i]) > 1 else ("SINGLE",)):
            for setName in (("devel", "test", "cafa") if (cafaTargets != "skip") else ("devel", "test")):
                combination = combinations[i][:]
                print
                print "***", "Evaluating", combination, "predictions for set '" + setName + "' using mode '" + mode + "'", "***"
                combinePred(task.proteins, combination, combKey, mode, limitToSets=[setName])
                examples = evaluateFile.makeExamples(task.proteins, limitTerms=limitTerms, limitToSets=[setName], predKey=combKey)
                loading.vectorizeExamples(examples, None, sparseLabels=True)
                try:
                    results = evaluation.evaluate(examples["labels"], examples["predictions"], examples, terms=None, averageOnly=True, noAUC=True)
                    print "Average for", str(combination) + "/" + setName + "/" + mode + ":", evaluation.metricsToString(results["average"])
                except TypeError as e:
                    print "WARNING! Cannot evaluate results."
                    print e
                if useOutFiles:
                    combString = "-".join(combination)
                    outPath = os.path.join(outDir, "-".join([combString, setName, mode, "ensemble"]) + ".tsv.gz")
                    evaluation.saveProteins(task.proteins, outPath, limitTerms=limitTerms, limitToSets=[setName], predKey=combKey) #pass#evaluation.saveResults(data, outStem, label_names, negatives)
                clearKeys(task.proteins, [combKey, combConfKey, combKey + "_sources"])

if __name__=="__main__":       
    from optparse import OptionParser
    optparser = OptionParser(description="Ensemble")
    optparser.add_option("-p", "--dataPath", default=os.path.expanduser("~/data/CAFA3/data"), help="Data directory")
    optparser.add_option("-i", "--inputs", default=None, help="Comma separated list of inputs in name:type:path format. Type is one of 'cls' or 'bl' (classifier or baseline). The first input must contain the gold labels.")    
    #optparser.add_option("-a", "--nnInput", default=None, help="Neural network predictions tsv.gz file")
    #optparser.add_option("-b", "--clsInput", default=None, help="Classifier predictions tsv.gz file")
    #optparser.add_option("-c", "--baselineInput", default=None, help="Baseline directory")
    optparser.add_option("-o", "--outDir", default=None, help="Output directory")
    optparser.add_option("-f", "--baselineCutoff", default=1, type=int, help="Cutoff for BLAST baseline predictions. Value in range 1-10, 1 for all values.")
    optparser.add_option("-m", "--modes", default="AND,OR", help="Input modes, comma-separated list of 'AND' and 'OR'")
    optparser.add_option("-t", "--terms", default=5000, type=int, help="The number of top most common GO terms to use as labels (will override task default)")
    optparser.add_option("-w", "--write", default=False, action="store_true", help="Write output files")
    optparser.add_option("--clear", default=False, action="store_true", help="Remove the output directory if it already exists")
    optparser.add_option("--targets", default="skip", help="How to include the CAFA target proteins, one of 'skip', 'overlap' or 'separate'. Default is 'skip'.")
    optparser.add_option("--task", default="cafa3")
    (options, args) = optparser.parse_args()
    
    combine(dataPath=options.dataPath, inputs=options.inputs, outDir=options.outDir,
            cafaTargets=options.targets, modes=options.modes, baselineCutoff=options.baselineCutoff,
            numTerms=options.terms, clear=options.clear, useOutFiles=options.write, taskName=options.task)