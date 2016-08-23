import sys, os
import argparse
from subprocess import check_call
import glob
import re
from itertools import izip
from mungo.fasta import FastaReader
from collections import defaultdict, Counter
import itertools

HMMERSEARCH = "hmmsearch"
DBLA_HMM = "/home/users/allstaff/tonkin-hill.g/DBLalpha/data/Rask_ref/atag.hmm"



def filterWithHMMER(inputfile, prefix, cpu, verbose):

  hmmer_cmd = (HMMERSEARCH
    + " -o /dev/null"
    + " --domT 80"
    + " --domtblout " + "hmmerDBLalphaSearch.txt"
    + " --cpu " + str(cpu)
    + " " + DBLA_HMM
    + " " + inputfile)

  if verbose:
    print "running... ", hmmer_cmd

  check_call(hmmer_cmd, shell=True)

  #Now run through and get DBLalpha seqs
  keep = set()
  with open("hmmerDBLalphaSearch.txt", 'rU') as infile:
    for line in infile:
      if line[0]=="#": continue
      keep.add(line.split()[0])

  with open(prefix + "_DBLa_cleaned.fasta", 'w') as dblfile:
    with open(prefix + "_NOT_dblalpha.fasta", 'w') as contamfile:
      for h,s in FastaReader(inputfile):
        if h in keep:
          dblfile.write(">" + h + "\n" + s + "\n")
        else:
          contamfile.write(">" + h + "\n" + s + "\n")


for f in glob.glob("/home/users/allstaff/tonkin-hill.g/DBLalpha/verificationOfData/controls/supp_material/demultiplexed_sequences/*.fasta"):
	filterWithHMMER(f, "/home/users/allstaff/tonkin-hill.g/DBLalpha/verificationOfData/controls/supp_material/hmmer_cleaned_sequeces/"+os.path.splitext(os.path.basename(f))[0], 30, True)
