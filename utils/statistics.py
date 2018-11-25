import os
import numpy
from collections import defaultdict
from sklearn.feature_extraction.dict_vectorizer import DictVectorizer

# From https://stackoverflow.com/questions/37604289/tkinter-tclerror-no-display-name-and-no-display-environment-variable
import matplotlib
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

def getCategories(sets):
    categories = sets
    if len(sets) == 1:
        categories += [sets[0] + " only"]
    return categories

def getCounts(examples):
    lengths = [len(examples[x]) for x in ("sets", "ids")] + [examples[x].shape[0] for x in ("features", "labels")]
    assert len(set(lengths)) == 1, lengths
    data = {}
    print "Counting labels"
    counts = defaultdict(list)
    for sets, labels in zip(examples["sets"], examples["labels"]):
        numLabels = numpy.count_nonzero(labels)
        #category = ",".join(sets)
        categories = getCategories(sets)
        for category in categories:
            counts[category].append(numLabels)
    data["labels"] = counts
    print "Counting features"
    dv = DictVectorizer(sparse=True)
    dv.feature_names_ = examples["feature_names"]
    decoded = dv.inverse_transform(examples["features"])
    for sets, features in zip(examples["sets"], decoded):
        categories = getCategories(sets)
        tags = [x.split(":")[0] for x in features.keys()]
        tagCounts = {x:tags.count(x) for x in set(tags)}
        for tag in tagCounts.keys():
            if tag not in data:
                data[tag] = defaultdict(list)
            counts = data[tag]
            for category in categories:
                counts[category].append(tagCounts[tag])
        
#         feature_indices = numpy.nonzero(features)
#         for i in feature_indices:
#             name = feature_names[i]
#             if ":" in name:
#                 tag = "coverage_" + name.split(":")[0]
#                 if tag not in data:
#                     data[tag] = defaultdict(list)
#                 counts = data[tag]
#                 for category in sets:
#                     counts[category] += 1
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
            ax.plot(x, values, label=category + " n=" + str(len(values)))
        plt.grid()
        plt.legend(fontsize = 'x-small')
        plt.ylabel("labels")
        plt.xlabel("proteins")
        if outDir != None:
            figPath = os.path.join(outDir, key + ".pdf")
            print "Saving figure to", figPath
            plt.savefig(figPath)
    
def makeStatistics(examples, outDir):
    data = getCounts(examples)
    makeGraphs(data, outDir)