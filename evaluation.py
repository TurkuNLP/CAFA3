import csv
import numpy as np
import gzip
from sklearn.metrics import roc_auc_score, f1_score, precision_score, recall_score

def evaluate(labels, predicted, examples, terms=None, averageOnly=False, average="micro"):
    print "Evaluating the predictions"
    results = {}
    print "Calculating average scores"
    results["average"] = {"id":"average", "ns":None, "name":None, "auc":0, "tp":None, "fp":None, "fn":None, "tn":None}
    try:
        results["average"]["auc"] = roc_auc_score(labels, predicted, average="micro")
    except ValueError as e:
        print e
    results["average"]["fscore"] = f1_score(labels, predicted, average=average)
    results["average"]["precision"] = precision_score(labels, predicted, average=average)
    results["average"]["recall"] = recall_score(labels, predicted, average=average)
    if averageOnly:
        return results
    
    print "Calculating label scores"
    label_names = examples["label_names"]
    label_size = examples.get("label_size")
    label_args = examples.get("label_args")
    try:
        aucs = roc_auc_score(labels, predicted, average=None)
    except ValueError as e:
        print e
        aucs = [0] * len(label_names)
    fscores = f1_score(labels, predicted, average=None)
    precisions = precision_score(labels, predicted, average=None)
    recalls = recall_score(labels, predicted, average=None)
    lengths = [len(x) for x in (aucs, fscores, precisions, recalls, label_names)]
    assert len(set(lengths)) == 1, lengths
    for auc, fscore, precision, recall, label_name in zip(aucs, fscores, precisions, recalls, label_names):
        assert label_name not in results
        result = {"id":label_name, "ns":None, "name":None, "auc":auc, "precision":precision, "recall":recall, "fscore":fscore, "tp":0, "fp":0, "fn":0, "tn":0}
        results[label_name] = result
        if label_size != None:
            result["label_size"] = label_size[label_name]
        if label_args != None:
            result["label_args"] = label_args[label_name]
        if terms != None and label_name in terms:
            term = terms[label_name]
            result["ns"] = term["ns"]
            result["name"] = term["name"]
    print "Counting label instances"
    stats = {x:{"tp":0, "fp":0, "fn":0, "tn":0} for x in label_names}
    label_indices = range(len(label_names))
    for gold, pred in zip(labels, predicted):
        for i in label_indices:
            stats[label_names[i]][getMatch(gold[i], pred[i])] += 1
    for key in stats:
        results[key].update(stats[key])
    return results

def resultIsBetter(original, new, key="average"):
    if new[key]["fscore"] != original[key]["fscore"]:
        return new[key]["fscore"] > original[key]["fscore"]
    else:
        return new[key]["auc"] > original[key]["auc"]
        
# def countInstances(labels, predicted):
#     tp = labels.multiply(predicted)
#     tp = labels.multiply(predicted)

def metricsToString(result, style="%.3f"):
    return "a/f|p/r|tp/fp/tn/fn = " + style % result["auc"] + "/" + style % result["fscore"] + "|" + style % result["precision"] + "/" + style % result["recall"] \
        + "|" + "/".join([str(result.get(x, "-")) for x in ("tp", "fp", "tn", "fn")])

def getResultsString(results, maxNumber=None, skipIds=None):
    count = 0
    s = ""
    for result in sorted(results.values(), key=lambda x: x["auc"], reverse=True):
        if skipIds != None and result["id"] in skipIds:
            continue
        s += metricsToString(result) + " " + str([result.get("id"), result.get("ns"), result.get("label_size"), result.get("name")]) + "\n"
        count += 1
        if count > maxNumber:
            break
    return s

def saveResults(data, outStem, label_names, negatives=False):
    print "Writing results to", outStem + "-results.tsv"
    with open(outStem + "-results.tsv", "wt") as f:
        dw = csv.DictWriter(f, ["auc", "fscore", "precision", "recall", "tp", "fp", "tn", "fn", "id", "label_size", "ns", "name", "label_args"], delimiter='\t')
        dw.writeheader()
        dw.writerow(data["results"]["average"])
        results = [x for x in data["results"].values() if x["id"] != "average"]
        dw.writerows(sorted(results, key=lambda x: x["auc"], reverse=True))
    savePredictions(data, label_names, outStem + "-predictions.tsv.gz", negatives=negatives)
    print "Writing ids to", outStem + "-ids.tsv"
    with open(outStem + "-ids.tsv", "wt") as f:
        dw = csv.DictWriter(f, ["id", "cafa_ids", "gold", "predicted"], delimiter='\t')
        dw.writeheader()
        dw.writerows([{"id":protId, "cafa_ids":",".join(cafa_ids), "gold":np.count_nonzero(gold), "predicted":np.count_nonzero(pred)} for protId, cafa_ids, gold, pred in zip(data["ids"], data["cafa_ids"], data["gold"], data["predicted"])]) 

def getMatch(gold, predicted):
    if gold == predicted:
        return "tp" if (gold == 1) else "tn"
    else:
        return "fn" if (gold == 1) else "fp"

def savePredictions(data, label_names, outPath, negatives=False):
    print "Writing predictions to", outPath
    keys = ["ids", "gold", "predicted", "cafa_ids"]
    hasProbabilities = data.get("probabilities") != None
    if hasProbabilities:
        lengths = [len(data["probabilities"]), len(label_names)]
        assert len(set(lengths)) == 1, lengths #keys += ["probabilities"]
    lengths = [len(data[x]) for x in keys]
    assert len(set(lengths)) == 1, lengths
    label_indices = range(len(label_names))
    
#     n_samples = data["probabilities"][0].shape[0]
#     label_names_array = np.array(label_names)
#     n_outputs = len(label_names_array)
#     predictions = np.zeros((n_samples, n_outputs))
#     for k in range(n_outputs):
#         predictions[:, k] = label_names_array[k].take(np.argmax(data["probabilities"][k], axis=1), axis=0)
    
    with gzip.open(outPath, "wt") as f:
        dw = csv.DictWriter(f, ["id", "label_index", "label", "predicted", "confidence", "gold", "match", "cafa_ids"], delimiter='\t')
        dw.writeheader()
        rows = []
        for i in range(len(data["ids"])):
            gold = data["gold"][i]
            pred = data["predicted"][i]
            cafa_ids = ",".join(data["cafa_ids"][i])
            for labelIndex in label_indices:
                goldValue = gold[labelIndex]
                predValue = int(pred[labelIndex])
                if negatives or (goldValue == 1 or predValue == 1):
                    row = {"id":data["ids"][i], "label_index":labelIndex, "label":label_names[labelIndex], "gold":goldValue, "predicted":predValue, "cafa_ids":cafa_ids}
                    row["match"] = getMatch(goldValue, predValue)
                    row["confidence"] = max(data["probabilities"][labelIndex][i]) if hasProbabilities else None #data["probabilities"][labelIndex][i] if hasProbabilities else None
                    #row["pred2"] = predictions[i][labelIndex]
                    rows.append(row)
            if len(rows) >= 100000:
                dw.writerows(rows)
                rows = []
        if len(rows) >= 0:
            dw.writerows(rows)