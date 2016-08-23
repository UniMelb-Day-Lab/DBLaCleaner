from mungo.fasta import FastaReader
import glob
import os, sys
from subprocess import check_call

USEARCH = "/home/users/allstaff/tonkin-hill.g/DBLalpha/third_party/usearch8.1.1831"
HMMERSEARCH = "hmmsearch"
DBLA_HMM = "/home/users/allstaff/tonkin-hill.g/DBLalpha/data/Rask_ref/atag.hmm"

def clusterWithUsearchConserved(combinedFile, per_id=0.96, min_conserved=3, verbose=True):
  #intially remove singletons
  usearch_cmd = (USEARCH
    + " -derep_prefix"
    + " " + combinedFile
    + " -fastaout " + combinedFile[:-6] + "_unique.fasta"
    + " -minuniquesize 2"
    + " -sizeout")
  if verbose:
    print "running... ", usearch_cmd
  check_call(usearch_cmd, shell=True)
  #cluster unique reads and annotate with size
  usearch_cmd = (USEARCH
    + " -cluster_fast"
    + " " + combinedFile[:-6] + "_unique.fasta"
    + " -centroids " + combinedFile[:-6] + "_centroids.fasta"
    + " -sort size"
    + " -id " + str(per_id))
  if verbose:
    print "running... ", usearch_cmd
  check_call(usearch_cmd, shell=True)
  #map all reads back to otus
  usearch_cmd = (USEARCH
    + " -usearch_global"
    + " " + combinedFile
    + " -db " + combinedFile[:-6] + "_centroids.fasta"
    + " -strand plus"
    + " -id " + str(per_id)
    + " -dbmatched " + combinedFile[:-6] + "_centroids_withSize.fasta"
    + " -sizeout"
    + " -otutabout " + combinedFile[:-6] + "_SingletonsFiltered_id" + str(per_id) + "_otutab.txt"
    )
  if verbose:
    print "running... ", usearch_cmd
  check_call(usearch_cmd, shell=True)
  return

for f in glob.glob("combined_*"):
  for pid in range(90,100,1):
    clusterWithUsearchConserved(f, pid/100.0)

