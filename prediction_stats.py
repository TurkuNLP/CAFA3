"""
Sanity check stats for test and cafa target predictions.
"""
import os
import sys
from collections import defaultdict
from scipy.stats import kendalltau, spearmanr

from evaluate import read_predictions

def stats(model_dir):
    test_pred = read_predictions(os.path.join(model_dir, 'test_pred.tsv.gz'))
    target_pred = read_predictions(os.path.join(model_dir, 'cafa_targets.tsv.gz'))

    print "Test set"
    test_gos = prediction_stats(test_pred)
    
    print "\nTarget set"
    target_gos = prediction_stats(target_pred)

    all_gos = list(set(test_gos.keys() + target_gos.keys()))

    # FIXME: I'm brain dead ATM, is this calculation correct?
    
    test_rank = [test_gos.get(go, 0) for go in all_gos]
    target_rank = [target_gos.get(go, 0) for go in all_gos]
    
    tau, pv = kendalltau(test_rank, target_rank)
    
    print "\nKendall's tau: %s, p-value: %s" % (tau, pv)
    
    print "\nSpearman's rho: %s, p-value: %s" % spearmanr(test_rank, target_rank)
    
    
def prediction_stats(predictions):
    proteins = len(predictions.keys())
    prediction_count= sum([len(v) for v in predictions.values()])
    gos = defaultdict(int)
    for prot_id, go_ids in predictions.items():
        for go in go_ids:
            gos[go] += 1
    
    unique_gos = len(gos.keys())
    
    print 'Unique proteins: %s' % proteins
    print 'Predictions in total: %s, average: %s' % (prediction_count, float(prediction_count)/proteins)
    print 'Unique GO concepts: %s' % unique_gos
    
    return gos
    
    

if __name__ == '__main__':
    stats(sys.argv[1])