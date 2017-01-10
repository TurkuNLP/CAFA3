import os
import numpy
from collections import defaultdict
import matplotlib.pyplot as plt

def getCounts(examples):
    lengths = [len(examples[x]) for x in ("sets", "ids")] + [examples[x].shape[0] for x in ("features", "labels")]
    assert len(set(lengths)) == 1, lengths
    counts = defaultdict(list)
    for sets, labels in zip(examples["sets"], examples["labels"]):
        example = {}
        numLabels = numpy.count_nonzero(labels)
        category = ",".join(sets)
        counts[category].append(numLabels)
    return counts

def makeGraph(counts, outPath):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for key in sorted(counts.keys()):
        values = sorted(counts[key])
        total = float(len(values))
        x = [x / total for x in reversed(range(len(values)))]
        ax.plot(x, values)
    plt.grid()
    plt.legend()
    plt.ylabel("labels")
    plt.xlabel("proteins")
    if outPath != None:
        plt.savefig(outPath)
    
def makeStatistics(examples, outDir):
    counts = getCounts(examples)
    makeGraph(counts, os.path.join(outDir, "labels.pdf"))