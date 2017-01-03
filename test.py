import sys, os
import gzip
from collections import defaultdict
import csv
from sklearn.feature_extraction import DictVectorizer
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.grid_search import GridSearchCV
from sklearn.ensemble.forest import ExtraTreesClassifier
import operator

def loadTerms(inPath, proteins):
    print "Loading terms from", inPath
    counts = defaultdict(int)
    with gzip.open(inPath, "rt") as f:
        tsv = csv.reader(f, delimiter='\t')
        for row in tsv:
            protId, goTerm, evCode = row
            protein = proteins[protId]
            if "terms" not in protein:
                protein["terms"] = {}
            protein["terms"][goTerm] = evCode
            counts[goTerm] += 1
    return counts

def loadSequences(inPath, proteins):
    print "Loading sequences from", inPath
    with gzip.open(inPath, "rt") as f:
        header = None
        for line in f:
            if header == None:
                assert line.startswith(">")
                header = line[1:].strip()
            else:
                proteins[header]["seq"] = line.strip()
                header = None
            #print seq.id, seq.seq

def buildExamples(proteins, limit=None, limitTerms=None):
    print "Building examples"
    y = []
    X = []
    mlb = MultiLabelBinarizer()
    dv = DictVectorizer(sparse=True)
    protIds = sorted(proteins.keys())
    if limit:
        protIds = protIds[0:limit]
    for protId in protIds:
        #print "Processing", protId
        protein = proteins[protId]
        # Build features
        features = {}
        seq = protein["seq"]
        for i in range(len(seq)-3):
            feature = seq[i:i+3]
            features[feature] = 1
        X.append(features)
        # Build labels
        labels = protein["terms"].keys()
        if limitTerms:
            labels = [x for x in labels if x in limitTerms]
        labels = sorted(labels)
        y.append(labels)
        #print features
    y = mlb.fit_transform(y)
    X = dv.fit_transform(X)
    return y, X

def getTopTerms(counts, num=1000):
    return sorted(counts.items(), key=operator.itemgetter(1), reverse=True)[0:num]

def classify(y, X, verbose=3, n_jobs = 1, scoring = "f1_micro"):
    clf = GridSearchCV(ExtraTreesClassifier(), {"n_estimators":[1,2,10,50,100]}, verbose=verbose, n_jobs=n_jobs, scoring=scoring)
    clf.fit(X, y)
    print "Best params", (clf.best_params_, clf.best_score_)

def run(dataPath):
    proteins = defaultdict(lambda: dict())
    loadSequences(os.path.join(options.dataPath, "Swiss_Prot", "Swissprot_sequence.tsv.gz"), proteins)
    counts = loadTerms(os.path.join(options.dataPath, "Swiss_Prot", "Swissprot_evidence.tsv.gz"), proteins)
    print "Proteins:", len(proteins)
    topTerms = getTopTerms(counts, 100)
    print "Most common terms:", topTerms
    print proteins["14310_ARATH"]
    y, X = buildExamples(proteins, 10000, set([x[0] for x in topTerms]))
    #print y
    #print X
    classify(y, X)

if __name__=="__main__":       
    from optparse import OptionParser
    optparser = OptionParser(description="")
    optparser.add_option("-p", "--dataPath", default=os.path.expanduser("~/data/CAFA3"), help="")
    (options, args) = optparser.parse_args()
    
    #proteins = de
    #importProteins(os.path.join(options.dataPath, "Swiss_Prot", "Swissprot_sequence.tsv.gz"))
    run(options.dataPath)