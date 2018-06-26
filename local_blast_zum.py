#  local_blast_zum.py
#  
#  Run on Python3

#  Created by Alice on 2018-06-26.
#

import argparse
from os import listdir
from os.path import isdir, isfile, join
from Bio.Blast.Applications import NcbiblastnCommandline

def main():
    args = get_arguments()
    output_file = "temp_out.txt"
    file_out = "res.txt"
    
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

    hits = {}
    for file in files:
        blastn_cmd=NcbiblastnCommandline(query=file, db="nematodeDB",
                                 max_target_seqs=1, gapopen=2, gapextend=3,
                                 outfmt="'6 qseqid sseqid stitle pident evalue length qstart qend mismatch gapopen gaps'", out=output_file)
        stdout, stderr = blastn_cmd()
        hits = summarize_blast_results(output_file, hits)

    save_2_file(hits,file_out)

def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="the input fasta files")
    parser.add_argument("-v", "--verbose", help="print more info",
                        action="store_true")
    return parser.parse_args()

def summarize_blast_results(results_file, hits):
    open(results_file) as res_file:
        for line in res_file:
            query_res = line.split('\t')
            hit_name = query_res[2]
            if hit_name in hits:
                hits[hit_name] += 1
            else:
                hits[hit_name] = 1
        res_file.close()

    return hits

def save_2_file(hits, file):
    tuple_hits = [ (organism, cnt) for organism, cnt in hits.items()]
    sorted_hits = sorted(tuple_hits, key=lambda x: x[1], reverse=True)
    sorted_hits_str = ["{}\t{}\n".format(organism, cnt) for organism, cnt in sorted_hits]

    filehandle = open(result_file, "w")
    filehandle.writelines(sorted_hits_str)
    filehandle.close()

main()
