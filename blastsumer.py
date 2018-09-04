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
from match_db import Match_db
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from numpy import log10

def main():
    args = get_arguments()
    result_file = "temp_out.txt"
    
    if isdir(args.input):
        files = [file for file in listdir(args.input)
                 if isfile(join(args.input, file)) and
                 file.split('.')[1]=="fa"]
    elif isfile(args.input):
        files = [args.input]
    else:
        raise NameError('Input file or directory does not exist')

    if args.verbose:
        print("\n---- Loaded input files ----")
        print(*files, sep='\n')

    hits ={}
    for file in files:
        if args.verbose:
            print("\nBlasting file: {}".format(file))
        blastn_cmd=NcbiblastnCommandline(query=file, db="../larvkult_1508/nematodeDB",
                                 max_target_seqs=1, gapopen=2, gapextend=3,
                                 outfmt="'6 qseqid sseqid stitle pident evalue length qstart qend mismatch gapopen gaps'", out=result_file)
        stdout, stderr = blastn_cmd()
        hits = summarize_blast_results(result_file, hits, args.pid, args.eval)

    save_2_file(hits,args.output, args.verbose)
    plot_results(hits, args.output)



def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="the input fasta files")
    parser.add_argument("-v", "--verbose", help="print more info",
                        action="store_true")
    parser.add_argument("-o", "--output", help="name of output directory",
                        required=True)
    parser.add_argument("--pid", help="Threshold of percentage identity of hits",
                        default=80.0, type=float)
    parser.add_argument("--eval", help="Threshold of e-val of hits",
                        default=1e-50, type=float)
    args = parser.parse_args()
    if args.output[len(args.output)-1]=='/':
        args.output = args.output[len(args.output)-2]
    if not exists(args.output):
        makedirs(args.output)
    return args

def summarize_blast_results(results_file, hits, perc_id_thresh, e_val_thresh):
    with open(results_file) as res_file:
        for line in res_file:
            query_res = line.split('\t')
            short_name = query_res[1]
            perc_id = float(query_res[3])
            e_val = float(query_res[4])
            if perc_id < perc_id_thresh or e_val > e_val_thresh:
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

    remove(results_file)
    return hits

def save_2_file(hits, output_dir, verbosity):
    output_file = output_dir + "/BLASTsumer_results.txt"
    tuple_hits = [ [match.name, match.count, mean(match.pid), mean(match.e_val),
                    mean(match.alg_len), mean(match.missmatch), mean(match.gaps),
                    mean(match.gap_openings)] for short_name, match in hits.items()]
    sorted_hits = sorted(tuple_hits, key=lambda x: x[1], reverse=True)
    sorted_hits_str = ["\t".join(list(map(str, hit)))+"\n" for hit in sorted_hits]
    
    if verbosity:
        print("\n---- Top 10 hits ----")
        print(*sorted_hits_str[:10], sep='')
        print("\nSaving results to: {}".format(output_dir))

    filehandle = open(output_dir, "w")
    filehandle.writelines(sorted_hits_str)
    filehandle.close()

def plot_results(hits, output_dir):
    sns.set(style="whitegrid")
    
    for name, hit in hits.items():
        if hit.count < 100:
            continue
        hit_dict = {'%id': hit.pid,
                    'length': hit.alg_len,
                    'missmatch': hit.missmatch,
                    'gaps': hit.gaps,
                    'gap openings': hit.gap_openings}
        hit_df = pd.DataFrame.from_dict(hit_dict)
        plot = sns.PairGrid(hit_df)
        plot = plot.map_upper(plt.scatter)
        plot = plot.map_lower(sns.kdeplot, cmap="Blues_d")
        plot = plot.map_diag(sns.distplot)
        fig = plot.fig
        fig.savefig(output_dir + '/' + name + "_result.png")
        print(hit_df.head())

main()
