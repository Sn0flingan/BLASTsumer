#  local_blast_zum.py
#  
#  Run on Python3

#  Created by Alice on 2018-06-26.
#

import argparse
from os import listdir, remove, makedirs
from os.path import exists, isdir, isfile, join
from statistics import mean
from Bio.Blast.Applications import NcbiblastnCommandline
import xml.etree.ElementTree as ET
from Bio.Blast import NCBIXML
from match_db import Match_db

def main():
    args = get_arguments()
    result_file = "temp_out.xml"
    
    if isdir(args.input):
        files = [file for file in listdir(args.input)
                 if isfile(join(args.input, file)) and
                 file.split('.')[1]=="fa"]
    elif isfile(args.input) and args.input.split(".")[-1]=="fa":
        files = [args.input]
    else:
        raise NameError('Input file or directory does not exist or is not .fa format')

    if args.verbose:
        print("\n---- Loaded input files ----")
        print(*files, sep='\n')

    result_file_2 = "temp_out2.txt"
    results = open(result_file_2, "w")
    hits ={}
    for file in files:
        if args.verbose:
            print("\nBlasting file: {}".format(file))
        blastn_cmd=NcbiblastnCommandline(query=file, db="../../larvkult_1508/nematodeDB",
                                 max_target_seqs=10, gapopen=2, gapextend=3,
                                 outfmt="5", out=result_file)
        stdout, stderr = blastn_cmd()

        #Parse results
        tree = ET.parse(result_file)
        root = tree.getroot()
        for read in root.iter('Iteration'):
            if int(read.find('Iteration_iter-num').text)%100==0:
                print("Read number: " + read.find('Iteration_iter-num').text)
            read_name = read.find('Iteration_query-def').text
            matches = {}
            for hit in read.iter('Hit'):
                hit_name = hit.find('Hit_def').text
                hit_len = float(hit.find('Hit_len').text)
                for hsp in hit.iter('Hsp'):
                    match_len = float(hsp.find('Hsp_identity').text)
                    perc_match = round(match_len/hit_len, 4)
                    identity = float(hsp.find('Hsp_identity').text)
                    e_val = hsp.find('Hsp_evalue').text
                    alg_len = hsp.find('Hsp_align-len').text
                    pid = (identity/float(alg_len))*100
                    gaps = int(hsp.find('Hsp_gaps').text)
                    missmatch = int(float(alg_len) - gaps - identity)
                    matches[hit_name] = [read_name, hit_name.split(" ")[0], hit_name,
                                         str(pid), e_val, str(alg_len), "0", "0", str(missmatch),
                                         "0", str(gaps), str(perc_match*100)]

            string_2_write = ""
            if len(matches)>0:
                best_hit = max(matches, key= lambda key: matches[key][-1])
                #print("{}% {}".format(matches[best_hit][-1], best_hit))
                string_2_write = matches[best_hit][0] + "\t" + matches[best_hit][1] \
                            + "\t" + matches[best_hit][2] + "\t" + matches[best_hit][3] \
                            + "\t" + matches[best_hit][4] + "\t" + matches[best_hit][5] \
                            + "\t" + matches[best_hit][6] +  "\t" + matches[best_hit][7] \
                            + "\t" + matches[best_hit][8] + "\t" + matches[best_hit][9] \
                            + "\t" + matches[best_hit][10] + "\t" + matches[best_hit][11] + "\n"
            else:
                string_2_write = read_name + "\tshortname\tlongname\t0.0\t1\t0\t0\t0\t0\t0\t0\t0"
                print("No hit found for read id: {}".format(read_name))
            results.write(string_2_write)


    results.close()
    hits = summarize_blast_results(result_file_2, hits, args.pid, args.eval)
    save_2_file(hits,args.output, args.verbose)


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="the input fasta files")
    parser.add_argument("-v", "--verbose", help="print more info",
                        action="store_true")
    parser.add_argument("-o", "--output", help="name of output file",
                        required=True)
    parser.add_argument("--pid", help="Threshold of percentage identity of hits",
                        default=60.0, type=float)
    parser.add_argument("--eval", help="Threshold of e-val of hits",
                        default=1e-10, type=float)
    args = parser.parse_args()
    if args.output.split(".")[-1]!='txt':
        args.output += ".txt"
    return args

def summarize_blast_results(results_file, hits, perc_id_thresh, e_val_thresh):
    with open(results_file) as res_file:
        for line in res_file:
            query_res = line.split('\t')
            short_name = query_res[1]
            perc_id = float(query_res[3])
            e_val = float(query_res[4])
            alg_len=int(query_res[5])
            if alg_len < 150 or perc_id < perc_id_thresh or e_val > e_val_thresh:
                short_name = 'None'
                query_res[2] = 'None'
            if short_name in hits:
                hits[short_name].add_read(read=query_res[0],
                                          pid=perc_id, alg_len=int(query_res[5]),
                                          e_val=e_val, missmatch=int(query_res[8]),
                                          gaps=int(query_res[10]),
                                          gaps_o=int(query_res[9]))
            else:
                hits[short_name] = Match_db(sn=short_name, name=query_res[2],
                                            pid=perc_id, alg_len=int(query_res[5]),
                                            e_val=e_val, missmatch=int(query_res[8]),
                                            gaps=int(query_res[10]),
                                            gaps_o=int(query_res[9]),
                                            read=query_res[0])
        res_file.close()

    #remove(results_file)
    return hits

def save_2_file(hits, output_file, verbosity):
    output_file = output_file
    tuple_hits = [ [match.name, match.count, mean(match.pid), mean(match.e_val),
                    mean(match.alg_len), mean(match.missmatch), mean(match.gaps),
                    mean(match.gap_openings)] for short_name, match in hits.items()]
    sorted_hits = sorted(tuple_hits, key=lambda x: x[1], reverse=True)
    sorted_hits_str = ["\t".join(list(map(str, hit)))+"\n" for hit in sorted_hits]
    
    if verbosity:
        print("\n---- Top 10 hits ----")
        print(*sorted_hits_str[:10], sep='')
        print("\nSaving results to: {}".format(output_dir))

    filehandle = open(output_file, "w")
    filehandle.writelines(sorted_hits_str)
    filehandle.close()

main()
