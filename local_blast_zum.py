#  local_blast_zum.py
#  
#  Run on Python3

#  Created by Alice on 2018-06-26.
#

import argparse

def main():
    args = get_arguments()

    if args.verbose:
        print("verbosity on")

    print(args.input)

def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="the input fasta files")
    parser.add_argument("-v", "--verbose", help="print more info",
                        action="store_true")
    return parser.parse_args()

main()
