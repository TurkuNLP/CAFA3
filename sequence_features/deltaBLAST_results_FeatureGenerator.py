#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os 
import sys
import gzip
from Bio.Blast import NCBIXML
import multiprocessing
import argparse


def func_os_GetAllFilesPathAndNameWithExtensionInFolder (FolderAddress, FileExtension, ProcessSubFoldersAlso = True):
    #IMPORTANT ... Extenstion should be like : "txt" , "a2"  ... WITHOUT DOT !    
    FILES = []; 
    if ProcessSubFoldersAlso:    
        for root, dirs, files in os.walk(FolderAddress):
            for file in files:
                if file.endswith("." + FileExtension):
                    FILES.append(os.path.join(root, file));
        return (FILES);
    else:         
        for file in os.listdir(FolderAddress):
            if file.endswith("." + FileExtension): #".txt" ;
                FILES.append(FolderAddress + file);
        return (FILES) ; 

    
def func_GetUniProtID_FromGI (hit_id , EVEXDBcon):
    #hit_id: gi|18202836|sp|Q9CQV8.3|1433B_MOUSE
    GI = hit_id.split ("|")[1]
    SQLSTR = "select symbol from map_gi2uniprot where gi = " + GI + " ;";
    cur = EVEXDBcon.cursor() 
    cur.execute (SQLSTR) 
    rows = cur.fetchall()
    cur.close()
  
    if len(rows) == 0:
        return set()  
    else:
        RES = set(); 
        for row in rows:
            RES.add (row[0])
        return RES 

def func_GetUniProtID_ACCESSION (hit_id , EVEXDBcon):
    #hit_id: gi|18202836|sp|Q9CQV8.3|1433B_MOUSE
    ACC = hit_id.split ("|")[3].split(".")[0] 
    SQLSTR = "select symbol from map_PID2uniprot where PID = '" + ACC + "' ;"
    cur = EVEXDBcon.cursor() 
    cur.execute (SQLSTR) 
    rows = cur.fetchall()
    cur.close()
  
    if len(rows) == 0:
        return set()  
    else:
        RES = set(); 
        for row in rows:
            RES.add (row[0])
        return RES 


def process_one_input_file (input_file , OUT_DIR):
    output_file = OUT_DIR + input_file.split("/")[-1].split(".xml.gz")[0] + ".features_tsv.gz" ;
    errlog_file = OUT_DIR + input_file.split("/")[-1].split(".xml.gz")[0] + ".errorlog.txt" ;

    print("processing       : " + input_file) 
    print("creating output  : " + output_file) 
    print("creating errorlog: " + errlog_file)

    inp_file_handle = gzip.open(input_file , 'rt') 
    out_file_handle = gzip.open(output_file, 'wt')
    log_file_handle = open(errlog_file, "wt")
    
    all_records = NCBIXML.parse(inp_file_handle) 
    cnt = 0 ; 
    for RECORD in all_records:
        for alignment in RECORD.alignments:
            for hsp in alignment.hsps:
                #cnt += 1 
                L = [] 
                #if cnt % 1000 == 0:
                #    print cnt ; 
                    
                            #Features : For each RECORD (Generally we have only 1 record)
                if len (RECORD.query.split ())>1:
                    # L.append (input_file.split("/")[-1].split(".")[0]);
                    L.append (RECORD.query.split ()[0]);#<BlastOutput_query-def>T96060004884 DPY30_HUMAN</BlastOutput_query-def> #UniprotID = DPY30_HUMAN
                else:
                    L.append (RECORD.query);
                    
                L.append (RECORD.query_id);#<BlastOutput_query-ID>90843</BlastOutput_query-ID>
                L.append (RECORD.query_length)#<Iteration_query-len>99</Iteration_query-len>
                L.append (RECORD.query_letters);#<BlastOutput_query-len>99</BlastOutput_query-len>
                                
                #Features : For each Alignment : EACH <Hit> ... and usually each <HIT> may have multiple <Hsp> ... we usually have 50 HSP
                # PARAM_UNIProtID_FromGI  = func_GetUniProtID_FromGI    (alignment.hit_id , EVEXDBcon) 
                # PARAM_UNIProtID_FromACC = func_GetUniProtID_ACCESSION (alignment.hit_id , EVEXDBcon) 
                # #hit_id: gi|18202836|sp|Q9CQV8.3|1433B_MOUSE
                PARAM_UNIProtID_FromXML = set() 
                tmp = alignment.hit_id.split("|") 
                if len (tmp) == 3:
                    PARAM_UNIProtID_FromXML.add (tmp[2])
                
                PARAM_UNIProtID = PARAM_UNIProtID_FromXML  
                if len(PARAM_UNIProtID) == 0:
                    ErrStr = RECORD.query_id + "\t" + alignment.hit_id + "\t" + "GI: " + alignment.hit_id.split ("|")[1] + "\n" ; 
                    log_file_handle.write (ErrStr) 
                    continue
                else:
                    PARAM_UNIProtID = ",".join (PARAM_UNIProtID)
                    
                L.append (PARAM_UNIProtID);# --> GI --> UniprotID 
                L.append (alignment.accession);#<Hit_accession>XP_005815176</Hit_accession>
                L.append (alignment.length);#<Hit_len>104</Hit_len>
                #L.append (alignment.hit_id);#<Hit_id>gi|551527403|ref|XP_005815176.1|</Hit_id>
                #L.append (alignment.hit_def);#<Hit_def>PREDICTED: protein dpy-30 homolog [Xiphophorus maculatus]</Hit_def>
    
                #Features : For each hsp : <hsp>
                L.append (hsp.align_length);#<Hsp_align-len>98</Hsp_align-len>
                L.append (hsp.bits) ;#<Hsp_bit-score>160.614</Hsp_bit-score>
                L.append (hsp.score);#<Hsp_score>405</Hsp_score>
                L.append (hsp.expect);# EVALUE : <Hsp_evalue>1.74162e-48</Hsp_evalue>
                L.append (hsp.query_start);#<Hsp_query-from>2</Hsp_query-from>
                L.append (hsp.query_end);#<Hsp_query-to>99</Hsp_query-to>
                L.append (hsp.sbjct_start);#<Hsp_hit-from>7</Hsp_hit-from>
                L.append (hsp.sbjct_end);#<Hsp_hit-to>104</Hsp_hit-to>
                L.append (hsp.frame[0]);#<Hsp_query-frame>0</Hsp_query-frame>
                L.append (hsp.frame[1]);#<Hsp_hit-frame>0</Hsp_hit-frame>
                L.append (hsp.identities);#<Hsp_identity>74</Hsp_identity>
                L.append (hsp.positives);#<Hsp_positive>92</Hsp_positive>
                L.append (hsp.gaps);#<Hsp_gaps>0</Hsp_gaps>
                
                out_file_handle.write ("\t".join(str(x) for x in L) + "\n")
                
    inp_file_handle.close()
    out_file_handle.close()
    log_file_handle.close() 

    
def worker_thread (input_file):
    args = argument_parser()
    OUT_DIR = args.out_folder
    print("PROCESSING FILE: " , input_file)
    process_one_input_file (input_file , OUT_DIR)
    print("PROCESSING FILE: " , input_file , " ... DONE.")

    
def argument_parser():
    parser = argparse.ArgumentParser(description="processXML output from interproscan")
    parser.add_argument("-i", "--in_folder", type=str, help="full path to input file -- xml.gz file from deltaBlast run")
    parser.add_argument("-o", "--out_folder", type=str, help="full path to out folder file from deltaBlast run")
    args = parser.parse_args()
    return args

    
if __name__== "__main__":
    args = argument_parser()
    INP_FILES = func_os_GetAllFilesPathAndNameWithExtensionInFolder(args.in_folder, "gz") 
    p = multiprocessing.Pool(8)
    p.map(worker_thread, INP_FILES)
    p.close() 
    p.join()     
