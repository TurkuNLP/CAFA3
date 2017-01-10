import os
import numpy
from collections import defaultdict
import matplotlib.pyplot as plt

def getCounts(examples):
    lengths = [len(examples[x]) for x in ("sets", "ids")] + [examples[x].shape[0] for x in ("features", "labels")]
    assert len(set(lengths)) == 1, lengths
    data = {}
    print "Counting labels"
    counts = defaultdict(list)
    for sets, labels in zip(examples["sets"], examples["labels"]):
        numLabels = numpy.count_nonzero(labels)
        #category = ",".join(sets)
        for category in sets:
            counts[category].append(numLabels)
    data["labels"] = counts
    print "Counting features"
    feature_names = examples["feature_names"]
    for sets, features in zip(examples["sets"], examples["features"]):
        feature_indices = numpy.nonzero(features)
        for i in feature_indices:
            name = feature_names[i]
            if ":" in name:
                tag = "coverage_" + name.split(":")[0]
                if tag not in data:
                    data[tag] = defaultdict(list)
                counts = data[tag]
                for category in sets:
                    counts[category] += 1
    return dict(data)

def makeGraphs(data, outDir):
    for key in sorted(data.keys()):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        counts = data[key]
        for category in sorted(counts.keys()):
            values = sorted(counts[category])
            total = float(len(values))
            x = [x / total for x in reversed(range(len(values)))]
            ax.plot(x, values, label=key)
        plt.grid()
        plt.legend()
        plt.ylabel("labels")
        plt.xlabel("proteins")
        if outDir != None:
            figPath = os.path.join(outDir, key + ".pdf")
            print "Saving figure to", figPath
            plt.savefig(figPath)
    
def makeStatistics(examples, outDir):
    data = getCounts(examples)
    makeGraphs(data, outDir)