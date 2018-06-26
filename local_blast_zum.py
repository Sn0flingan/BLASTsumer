#  local_blast_zum.py
#  
#  Run on Python3

#  Created by Alice on 2018-06-26.
#

import argparse
from os import listdir
from os.path import isdir, isfile, join

def main():
    args = get_arguments()
    
    if isdir(args.input):
        files = [file for file in listdir(args.input)
                 if isfile(join(args.input, file)) and
                 file.split('.')[1]=="fa"]
    elif isfile(args.input):
        files = [args.input]
    else:
        raise NameError('Input file or directory does not exist')
    
    print(files)

def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="the input fasta files")
    parser.add_argument("-v", "--verbose", help="print more info",
                        action="store_true")
    return parser.parse_args()

main()
