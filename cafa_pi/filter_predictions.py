"""
Removes motility predictions from Fungi. Only Bacteria should have motility annotations.
"""

import sys

import gzip
import networkx as nx

    
parent_file = gzip.open('./data/go_parent_term2_term.tsv.gz')
parent_data = parent_file.readlines()[1:] # First line has the column names

print('Building GO graph')
go_graph = nx.DiGraph()

for line in parent_data:
    parent, child = line.strip().split('\t')
    go_graph.add_edge(parent, child)

in_path = sys.argv[1]
out_path = sys.argv[2]

f = gzip.open(in_path, 'rt')
out = []
removed = 0
for line in f:
    if line.startswith('C'):
        go = line.split()[2]
        if not 'GO:0001539' in nx.ancestors(go_graph, go): # This is the root for the motility subtree in GO
            out.append(line)
        else:
            removed += 1
    else:
        out.append(line)

out_f = gzip.open(out_path, 'wt')
out_f.writelines(out)
print 'Removed %s annotations' % removed