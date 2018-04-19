"""
Converts our internal format to the oofficial CAFA PI format.
The output should have one file per organism per GO term, i.e. 4 files in total.
"""
import os
import sys
import gzip
import csv

GROUP_NAME = 'TurkuNLP'

input_path = sys.argv[1]
output_directory = sys.argv[2]
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

model_number = sys.argv[3]

KEYWORDS = ['convolutional neural network', 'random forest']

in_file = gzip.open(input_path, 'rt')
reader = csv.reader(in_file, delimiter='\t', quotechar='|')
reader.next() # Skip header line

predictions = {208963: {'GO:0001539': [], 'GO:0042710': []}, # Bacteria dude
               237561: {'GO:0042710': []}} # Fungi dude

for prot_id, label_index, label, predicted, confidence, cafa_id in reader:
    if prot_id.startswith('C'):
        predictions[237561][label].append((cafa_id, confidence))
    else:
        try:
            predictions[208963][label].append((cafa_id, confidence))
        except:
            print label
            import pdb; pdb.set_trace()

for organism in predictions:
    for go_predictions in predictions[organism]: # All proteins for a specific GO and specific organism
        out_file = open(os.path.join(output_directory, '%s_%s_%s_%s.txt' % (GROUP_NAME, model_number, organism, go.replace('GO:', ''))))
        out_file.write('AUTHOR\t%s\n' % GROUP_NAME)
        out_file.write('MODEL\t%s\n' % model_number)
        out_file.write('KEYWORDS\t%s\n' % ', '.join(KEYWORDS))
        
        for cafa_id, confidence in go_predictions:
            out_file.write('%s\t%s\n' % (cafa_id, confidence))
        
        out_file.write('END')