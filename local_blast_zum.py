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

    return

    for file in files:
        blastn_cmd=NcbiblastnCommandline(query=file, db="nematodeDB",
                                 max_target_seqs=1, gapopen=2, gapextend=3,
                                 outfmt="'6 qseqid sseqid stitle pident evalue length qstart qend mismatch gapopen gaps'", out=output_file)
        stdout, stderr = blastn_cmd()

def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="the input fasta files")
    parser.add_argument("-v", "--verbose", help="print more info",
                        action="store_true")
    return parser.parse_args()


main()
