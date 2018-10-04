#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse
import glob
import gzip

for filename in glob.glob('/home/sukaew/CAFA3/interproscan_result_features/*'):
    with gzip.open(filename, 'rb') as f:
        all_lines = f.readlines()
    with gzip.open(filename.replace('CAFA3/interproscan_result_features/', 'CAFA_PI/targetFiles/interproscan_features/'), 'rb') as f:
        new_lines = f.read()
    if [new_lines.split('\n')[0]] == ['']:
        new_text = all_lines[0].strip('\n') + ''.join(new_lines)
    else:
        new_text = ''.join(new_lines)
    with gzip.open(filename.replace('CAFA3/interproscan_result_features/', 'CAFA_PI/targetFiles/temp_interproscan_features/'), 'wb') as f:
       f.write(new_text)
    
