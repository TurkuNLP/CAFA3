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
        for p in parent_annotations:
            new_data.add((prot_id, p))
            
    prop_file = gzip.open('./data/Swissprot_propagated.tsv.gz', 'w')
    
    print 'Writing to file'
    for p_id, go_id in new_data:
        prop_file.write('%s\t%s\tNONE\n' % (p_id, go_id)) # FIXME: Evidence codes are lost for now

if __name__ == '__main__':
    propagate()