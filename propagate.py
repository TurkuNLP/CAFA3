"""
Propagates parent annotations for all sequences.
"""
import gzip
import networkx as nx

def propagate():
    ann_file = gzip.open('./data/Swissprot_evidence.tsv.gz')
    ann_data = ann_file.readlines()
    
    parent_file = gzip.open('./data/parent_term2term.txt.gz')
    parent_data = parent_file.readlines()[1:] # First line has the column names
    
    print 'Building GO graph'
    go_graph = nx.DiGraph()
    
    for line in parent_data:
        parent, child = line.strip().split('\t')
        go_graph.add_edge(parent, child)
    
    print 'Propagating annotations'
    new_data = set()
    for line in ann_data:
        prot_id, annotation, evidence = line.strip().split('\t')
        new_data.add((prot_id, annotation))
        try:
            parent_annotations = nx.ancestors(go_graph, annotation)
        except:
            print 'ERROR: %s not in GO graph (probably obsolete GO term)' % annotation
            continue
        for p in parent_annotations:
            new_data.add((prot_id, p))
            
    prop_file = gzip.open('./data/Swissprot_propagated.tsv.gz', 'w')
    
    print 'Writing to file'
    for p_id, go_id in new_data:
        prop_file.write('%s\t%s\tNONE\n' % (p_id, go_id)) # FIXME: Evidence codes are lost for now

def propagate_predictions(predictions):
    parent_file = gzip.open('./data/parent_term2term.txt.gz')
    parent_data = parent_file.readlines()[1:]
    go_graph = nx.DiGraph()
    
    for line in parent_data:
        parent, child = line.strip().split('\t')
        go_graph.add_edge(parent, child)
    
    parent_dict = {}
    
    new_pred = []
    not_found = 1
    for i, preds in enumerate(predictions):
        if i % 1000 == 0:
            print i
        pred_set = set()
        for pe in preds:
            pred_set.add(pe)
            try:
                if pe in parent_dict:
                    parent_annotations = parent_dict[pe]
                else:
                    parent_annotations = nx.ancestors(go_graph, pe)
                    parent_dict[pe] = parent_annotations
            except:
                #print 'PROPAGATION ERROR: %s not in GO graph (probably obsolete GO term)' % pe
                #import pdb; pdb.set_trace()
                not_found += 1
                continue
            pred_set.update(parent_annotations)
        new_pred.append(list(pred_set))
    print '%s labels not propagated' % not_found
    return new_pred

if __name__ == '__main__':
    propagate()