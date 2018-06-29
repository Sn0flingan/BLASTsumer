# BLAST sum:er
Simple python script to run BLAST on a local database over many sequences and sum up the results. Also includes proper output to analyse how different statistics such as % identity, alignment length, missmatches...etc. covariate.

To run the script:
```
python3 blastsumer.py -i <input file or directory in .fa (fasta) format> -o <output directory (does not need to exist)>
```

Optional parameters include:
- `-v` or `--verbose` add this argument for more verbose output in the terminal while running
- `--pid` this is a threshold for minimum % ID of alignment of BLAST hits. Default set to 80.0 (80%)
- `--eval` this is a threshold for maximum e-value of BLAST hits. Default set to 1e-50. 
