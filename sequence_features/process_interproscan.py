#!/usr/bin/env python
# -*- coding: utf-8 -*-


from xml.etree import cElementTree as ET
import argparse
import gzip
import glob
import re


def load_mapping(GO_map_folder):
    db_source = {}
    db_mapping = []
    for in_file in glob.glob(GO_map_folder + '*'):
        with gzip.open(in_file, "rt") as f:
            for line in f:
                if line[0] != "!" and "COG" not in line.strip("\n").split(" ; ")[-1]:
                    db_mapping.append([line.split(":")[1].split(" > ")[0].split(" ")[0], line.strip("\n").split(" ; ")[-1]])
                
            # db_mapping += [[line.split(":")[1].split(" > ")[0].split(" ")[0], line.strip("\n").split(" ; ")[-1]] for line in f if line[0] != "!" and "COG" not in line.strip("\n").split(" ; ")[-1]]
            db_source.setdefault(line.split(":")[1].split(" > ")[0].split(" ")[0], in_file.split("/")[-1].strip("2go"))
    db_dict = dict()
    for item in db_mapping:
        db_dict.setdefault(item[0], []).append(item[1])
    return db_dict, db_source

                
def extract_result_GO(protein_id, panther, name, db_dict, db_source):
    all_match = []
    header = ""
    for match in panther:
        detail = []
        db_id = []
        entry = []
        count = 0
        evalue = []
        detail_col = []
        evalue_col = []
        for elem in match.iter():
            elem_tag = elem.tag.split("}")[1]
            if elem_tag in ["{}-location".format(name)]:
                count += 1
                detail += [v for k, v in elem.attrib.items() if k not in ["familyName", "level"]]
                detail_col = [k for k, v in elem.attrib.items() if k not in ["familyName", "level"]]
            elif elem_tag in ["{}-match".format(name)]:
                evalue += [v for k, v in elem.attrib.items() if k not in ["familyName", "level"]]
                evalue_col = [k for k, v in elem.attrib.items() if k not in ["familyName", "level"]]
            elif elem_tag in ["go-xref"]:
                # print [v for k, v in elem.attrib.items() if k in ["db", "id"]]
                db_id.append([v for k, v in elem.attrib.items() if k in ["db", "id"]])
            elif elem_tag in ["signature", "entry", "model"]:
                try:
                    for go_match in db_dict[elem.attrib['ac']]:
                        db_id.append([db_source[elem.attrib['ac']], go_match])
                        # db_id.append([db_dict[elem.attrib['ac']], go_match])
                except:
                    pass
            elif elem_tag in ["pathway-xref"]:
                try:
                    for go_match in db_dict[elem.attrib['id']]:
                        db_id.append([db_source[elem.attrib['id']], go_match])
                        # db_id.append([db_dict[elem.attrib['ac']], go_match])
                except:
                    pass
        header = "\t".join(["protein_id", "db"] + evalue_col + detail_col)
        l_index = int(len(detail)/count)
        for item in db_id:
            i = 0
            while i < len(detail):
                # print(i, i+l_index, type(i), type(i+l_index))
                text = "\t".join([protein_id] + item + evalue + detail[i:i+l_index])
                if "GO" in text:
                    all_match.append(text)
                i += l_index                
    return header, "\n".join(all_match)

                        
def process_xml_GO(in_folder, out_folder, GO_map_folder):
    """
    """
    db_dict, db_source = load_mapping(GO_map_folder)
    all_analysis =  ["hmmer3", "panther", "patternscan", "phobius", "profilescan", "rpsblast", "superfamilyhmmer3", "coils", "hmmer2", "fingerprints", "blastprodom", "tmhmm"]
    library = []
    totalString = ""
    basedTag = "{http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5}"
    allDict = dict()
    header_dict = {"panther" : ["protein_id", "db", "GOid", "evalue", "score", "start", "end"],
                   "hmmer3": ["protein_id", "db", "GOid", "evalue", "score", "end", "env-end", "hmm-lengthe", "value", "start", "score", "hmm-end", "hmm-start", "env-start"],                   
                   "patternscan": ["protein_id", "db", "GOid", "start", "end"],                   
                   "phobius": ["protein_id", "db", "GOid", "start", "end"],                   
                   "profilescan": ["protein_id", "db", "GOid", "start", "score", "end"],                   
                   "rpsblast": ["protein_id", "db", "GOid", "start", "evalue", "score", "end"],
                   "superfamilyhmmer3": ["protein_id", "db", "GOid", "evalue", "start", "end"],
                   "coils": ["protein_id", "db", "GOid", "start", "end"],                   
                   "hmmer2": ["protein_id", "db", "GOid", "evalue", "score", "end", "hmm-length", "evalue", "start", "score", "hmm-end", "hmm-start"],
                   "fingerprints": ["protein_id", "db", "GOid", "evalue", "graphscan", "start", "score", "end", "motifNumber", "pvalue"],
                   "blastprodom": ["protein_id", "db", "GOid", "start", "score", "end"],
                   "tmhmm": ["protein_id", "db", "GOid", "start", "end"]}
    allXML = glob.glob(in_folder + "*.xml.gz")
    allXML.sort()
    for in_xml in allXML:
        outputFile = in_xml.split(".")[0] + ".txt"
        tree = ET.parse(gzip.open(in_xml))
        doc = tree.getroot() # get root element
        proteins = doc.findall(".//{}protein".format(basedTag))
        for protein in proteins:
            protein_id = protein.findall(".//{}xref".format(basedTag))[0].attrib['id']
            for analysis in all_analysis:
                header, result = extract_result_GO(protein_id, protein.findall(".//{}{}-match".format(basedTag, analysis)), analysis, db_dict, db_source)
                allDict.setdefault(analysis, []).append(result)
                # header_dict.setdefault(analysis, [header])
    for k, v in allDict.items():
        with gzip.open(out_folder + k + "_GO.tsv.gz", "wt") as x:
            split_v = [item for item in v if "" not in item.split("\t")]
            # print([item.split('\t') for item in v])
            # if len(split_v) > 0:
            #     len(split_v[0].split('\t')) == len(header[k])
            x.write('\t'.join(header_dict[k]) + '\n' + '\n'.join(split_v) + '\n')

            
def extract_result_noGO(protein_id, panther, name, db_dict, db_source):
    all_match = []
    header = ""
    for match in panther:
        detail = []
        db_id = []
        entry = []
        count = 0
        evalue = []
        detail_col = []
        evalue_col = []
        for elem in match.iter():
            elem_tag = elem.tag.split("}")[1]
            if elem_tag in ["{}-location".format(name)]:
                count += 1
                detail += [v for k, v in elem.attrib.items() if k not in ["familyName", "level"]]
                detail_col = [k for k, v in elem.attrib.items() if k not in ["familyName", "level"]]
            elif elem_tag in ["{}-match".format(name)]:
                evalue += [v for k, v in elem.attrib.items() if k not in ["familyName", "level"]]
                evalue_col = [k for k, v in elem.attrib.items() if k not in ["familyName", "level"]]
            elif elem_tag in ["signature", "entry"]:                
                db_id.append([v for k, v in elem.attrib.items() if k in ["ac"]])
            header = "\t".join(["protein_id", "db"] + evalue_col + detail_col)
        l_index = int(len(detail)/count)
        for item in db_id:
            i = 0
            while i < len(detail):
                text = "\t".join([protein_id] + item + evalue + detail[i:i+l_index])
                all_match.append(text)
                i += l_index
    return header, "\n".join(all_match)

                        
def process_xml_noGO(in_folder, out_folder, GO_map_folder):
    """
    """
    db_dict, db_source = load_mapping(GO_map_folder)
    all_analysis =  ["panther", "hmmer3", "patternscan", "phobius", "profilescan", "rpsblast", "superfamilyhmmer3", "coils", "hmmer2", "fingerprints", "blastprodom", "tmhmm"]
    library = []
    totalString = ""
    basedTag = "{http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5}"
    allDict = {}
    header_dict = {"panther": ["protein_id", "ac", "evalue", "score", "start", "end"],
                   "hmmer3": ["protein_id", "db", "evalue", "score", "end", "env-end", "hmm-lengthe", "value", "start", "score", "hmm-end", "hmm-start", "env-start"],
                   "patternscan": ["protein_id", "ac", "start", "end"],
                   "phobius": ["protein_id", "ac", "start", "end"],
                   "profilescan": ["protein_id", "ac", "start", "score", "end"],
                   "rpsblast": ["protein_id", "ac", "start", "evalue", "score", "end"],
                   "superfamilyhmmer3": ["protein_id", "db", "evalue", "start", "end"],
                   "coils": ["protein_id", "db", "start", "end"],                   
                   "hmmer2": ["protein_id", "ac", "evalue", "score", "end", "hmm-length", "evalue", "start", "score", "hmm-end", "hmm-start"],
                   "fingerprints": ["protein_id", "ac", "evalue", "graphscan", "start", "score", "end", "motifNumber", "pvalue"],
                   "blastprodom": ["protein_id", "ac", "start", "score", "end"],
                   "tmhmm": ["protein_id", "db", "start", "end"]}
    allXML = glob.glob(in_folder + "*.xml.gz")
    allXML.sort()
    for in_xml in allXML:
        outputFile = in_xml.split(".")[0] + ".txt"
        tree = ET.parse(gzip.open(in_xml))
        doc = tree.getroot() # get root element        
        proteins = doc.findall(".//{}protein".format(basedTag))
        for protein in proteins:
            protein_id = protein.findall(".//{}xref".format(basedTag))[0].attrib['id']            
            for analysis in all_analysis:
                header, result = extract_result_noGO(protein_id, protein.findall(".//{}{}-match".format(basedTag, analysis)), analysis, db_dict, db_source)
                allDict.setdefault(analysis, []).append(result)                
                # header_dict.setdefault(analysis, [header])
    for k, v in allDict.items():
        with gzip.open(out_folder + k + "_noGO.tsv.gz", "wt") as x:
            split_v = [item for item in v if "" not in item.split("\t")]
            x.write('\t'.join(header_dict[k]) + '\n' + '\n'.join(split_v) + '\n')
            
            
def check_header(out_folder):
    """
    check if the rows and header contain the same number of columns
    """
    for filename in glob.glob(out_folder + '*.gz'):
        with gzip.open(filename, 'rt') as f:
            for i, line in enumerate(f):
                if i == 0:
                    col_count = len(line.split('\t'))
                else:
                    if line != '\n':
                        assert len(line.split('\t')) == col_count #rows and header is not the same number of columns

                        
def argument_parser():
    parser = argparse.ArgumentParser(description="process XML output from interproscan")
    parser.add_argument("-i", "--input_folder", type=str, help="full path to input file -- xml.gz file from interproscan run")
    parser.add_argument("-o", "--out_folder", type=str, help="full path to out folder file from interproscan run")
    parser.add_argument("-m", "--GO_map_folder", type=str, help="full path to out folder file from interproscan run")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = argument_parser()    
    process_xml_GO(args.input_folder, args.out_folder, args.GO_map_folder)
    check_header(args.out_folder)
    process_xml_noGO(args.input_folder, args.out_folder, args.GO_map_folder)
    check_header(args.out_folder)
