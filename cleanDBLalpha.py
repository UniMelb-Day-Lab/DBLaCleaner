import sys, os
import argparse
from subprocess import check_call
import glob
import re
from itertools import izip
from mungo.fasta import FastaReader

PEAR = "pear"
FLEXBAR = "/home/users/allstaff/tonkin-hill.g/DBLalpha/third_party/flexbar_v2.5_linux64/flexbar-2.5.0/flexbar"
USEARCH = "/home/users/allstaff/tonkin-hill.g/DBLalpha/third_party/usearch8.1.1831"
HMMERSEARCH = "hmmsearch"
DBLA_HMM = "/home/users/allstaff/tonkin-hill.g/DBLalpha/data/Rask_ref/atag.hmm"

FORWARD_PRIMER = "GNNNGNAGTTTNGC"
REVERSE_PRIMER = "NNCCANTCNTNNNACCA"
# MID1-153
MIDIC = {"1":"ACGAGTGCGT","40":"TACGCTGTCT","2":"ACGCTCGACA","41":"TAGTGTAGAT","3":"AGACGCACTC","42":"TCGATCACGT","4":"AGCACTGTAG","43":"TCGCACTAGT","5":"ATCAGACACG","44":"TCTAGCGACT","6":"ATATCGCGAG","45":"TCTATACTAT","7":"CGTGTCTCTA","46":"TGACGTATGT","8":"CTCGCGTGTC","47":"TGTGAGTAGT","9":"TAGTATCAGC","10":"TCTCTATGCG","48":"ACAGTATATA","11":"TGATACGTCT","49":"ACGCGATCGA","12":"TACTGAGCTA","13":"CATAGTAGTG","50":"ACTAGCAGTA","14":"CGAGAGATAC","51":"AGCTCACGTA","15":"ATACGACGTA","52":"AGTATACATA","16":"TCACGTACTA","53":"AGTCGAGAGA","17":"CGTCTAGTAC","54":"AGTGCTACGA","18":"TCTACGTAGC","55":"CGATCGTATA","19":"TGTACTACTC","56":"CGCAGTACGA","20":"ACGACTACAG","57":"CGCGTATACA","21":"CGTAGACTAG","58":"CGTACAGTCA","22":"TACGAGTATG","59":"CGTACTCAGA","23":"TACTCTCGTG","60":"CTACGCTCTA","24":"TAGAGACGAG","61":"CTATAGCGTA","25":"TCGTCGCTCG","62":"TACGTCATCA","26":"ACATACGCGT","63":"TAGTCGCATA","27":"ACGCGAGTAT","64":"TATATATACA","28":"ACTACTATGT","65":"TATGCTAGTA","29":"ACTGTACAGT","66":"TCACGCGAGA","30":"AGACTATACT","67":"TCGATAGTGA","31":"AGCGTCGTCT","68":"TCGCTGCGTA","32":"AGTACGCTAT","69":"TCTGACGTCA","33":"ATAGAGTACT","70":"TGAGTCAGTA","34":"CACGCTACGT","71":"TGTAGTGTGA","35":"CAGTAGACGT","72":"TGTCACACGA","36":"CGACGTGACT","73":"TGTCGTCGCA","37":"TACACACACT","74":"ACACATACGC","38":"TACACGTGAT","75":"ACAGTCGTGC","39":"TACAGATCGT","76":"ACATGACGAC","77":"ACGACAGCTC","78":"ACGTCTCATC","79":"ACTCATCTAC","80":"ACTCGCGCAC","81":"AGAGCGTCAC","82":"AGCGACTAGC","83":"AGTAGTGATC","84":"AGTGACACAC","85":"AGTGTATGTC","86":"ATAGATAGAC","87":"ATATAGTCGC","88":"ATCTACTGAC","89":"CACGTAGATC","90":"CACGTGTCGC","91":"CATACTCTAC","92":"CGACACTATC","93":"CGAGACGCGC","94":"CGTATGCGAC","95":"CGTCGATCTC","96":"CTACGACTGC","97":"CTAGTCACTC","98":"CTCTACGCTC","99":"CTGTACATAC","100":"TAGACTGCAC","101":"TAGCGCGCGC","102":"TAGCTCTATC","103":"TATAGACATC","104":"TATGATACGC","105":"TCACTCATAC","106":"TCATCGAGTC","107":"TCGAGCTCTC","108":"TCGCAGACAC","109":"TCTGTCTCGC","110":"TGAGTGACGC","111":"TGATGTGTAC","112":"TGCTATAGAC","113":"TGCTCGCTAC","114":"ACGTGCAGCG","115":"ACTCACAGAG","116":"AGACTCAGCG","117":"AGAGAGTGTG","118":"AGCTATCGCG","119":"AGTCTGACTG","120":"AGTGAGCTCG","121":"ATAGCTCTCG","122":"ATCACGTGCG","123":"ATCGTAGCAG","124":"ATCGTCTGTG","125":"ATGTACGATG","126":"ATGTGTCTAG","127":"CACACGATAG","128":"CACTCGCACG","129":"CAGACGTCTG","130":"CAGTACTGCG","131":"CGACAGCGAG","132":"CGATCTGTCG","133":"CGCGTGCTAG","134":"CGCTCGAGTG","135":"CGTGATGACG","136":"CTATGTACAG","137":"CTCGATATAG","138":"CTCGCACGCG","139":"CTGCGTCACG","140":"CTGTGCGTCG","141":"TAGCATACTG","142":"TATACATGTG","143":"TATCACTCAG","144":"TATCTGATAG","145":"TCGTGACATG","146":"TCTGATCGAG","147":"TGACATCTCG","148":"TGAGCTAGAG","149":"TGATAGAGCG","150":"TGCGTGTGCG","151":"TGCTAGTCAG","152":"TGTATCACAG","153":"TGTGCGCGTG"}


class writeable_dir(argparse.Action):
  def __call__(self,parser, namespace, values, option_string=None):
      prospective_dir=values
      if not os.path.isdir(prospective_dir):
          raise argparse.ArgumentTypeError("writeable_dir:{0} is not a valid path".format(prospective_dir))
      if os.access(prospective_dir, os.W_OK):
          setattr(namespace,self.dest,prospective_dir)
      else:
          raise argparse.ArgumentTypeError("writeable_dir:{0} is not a writeable dir".format(prospective_dir))

class readable_file(argparse.Action):
  def __call__(self,parser, namespace, values, option_string=None):
      prospective_file=values
      if not os.path.isfile(prospective_file):
          raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid file path".format(prospective_file))
      if os.access(prospective_file, os.R_OK):
          setattr(namespace,self.dest,prospective_file)
      else:
          raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable file".format(prospective_file))

class MyError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def isValidBarcode(barcode):
  for l in barcode:
    if l not in ["A","C","G","T","N"]:
      return False
  return True

def catFiles(filenames, outputfile):
  with open(outputfile, 'w') as outfile:
      for fname in filenames:
          with open(fname, 'r') as infile:
              for line in infile:
                  outfile.write(line)

def checkFastqPaired(read1file, read2file, verbose):
  count = 0
  with open(read1file, 'r') as read1:
    with open(read2file, 'r') as read2:
      for r1, r2 in izip(read1, read2):
        if r1[0]=="@" or r2[0]=="@":
          count += 1
          if r1.split()[0]!=r2.split()[0]:
            MyError("Read file mismatch")

  if count<1:
    MyError("Empty read files!")

  if verbose:
    print "Checked ", count, " read pairs."

  return count


def convertDescToFlexBarcodes(descfile, outputfile, verbose):
  pairedMIDs = {}
  barcodes = set()

  with open(descfile, 'r') as infile:
    with open(outputfile, 'w') as outfile:
      for line in infile:

        #Check barcodes
        if line[0]=="#": continue
        line = line.strip().split()
        if line[1] in MIDIC:
          b1 = MIDIC[line[1]]
        elif isValidBarcode(line[1]):
          b1 = line[1]
        else:
          MyError("Invalid barcode in desc file")
        if line[2] in MIDIC:
          b2 = MIDIC[line[2]]
        elif isValidBarcode(line[2]):
          b2 = line[2]
        else:
          MyError("Invalid barcode in desc file")

        if ((line[1], line[2]) in pairedMIDs) or ((line[2], line[1]) in pairedMIDs):
          MyError("Duplicated MID id pair.")

        pairedMIDs[(line[1], line[2])] = line[0]

        #write barcode sequences with primers
        if b1 not in barcodes:
          outfile.write(">B_" + line[1] + "_forward\n")
          outfile.write(b1 + FORWARD_PRIMER +"\n")
          outfile.write(">B_" + line[1] + "_reverse\n")
          outfile.write(b1 + REVERSE_PRIMER +"\n")
          barcodes.add(b1)

        if b2 not in barcodes:
          outfile.write(">B_" + line[2] + "_forward\n")
          outfile.write(b2 + FORWARD_PRIMER +"\n")
          outfile.write(">B_" + line[2] + "_reverse\n")
          outfile.write(b2 + REVERSE_PRIMER +"\n")
          barcodes.add(b2)

  return outputfile, pairedMIDs

def demultiplexTrim(barcodefile, read1, read2, cpu, verbose):

  #Check read inputs
  if verbose:
    print "Checking paired end reads match..."

  total_paired_count = checkFastqPaired(read1, read2, verbose)

  opts = (" --barcode-threshold 1"
    + " --max-uncalled 15"
    + " --removal-tags"
    + " --barcode-trim-end LEFT"
    + " --min-read-length 100"
    # + " --pre-trim-phred 28"
    # + " --barcode-unassigned"
    + " --format sanger"
    + " --threads " + str(cpu)
    + " --barcodes " + barcodefile)

  #Run flexbar splitting on read1
  flexbar_cmd = (FLEXBAR
    + " --reads " + read1
    + " --reads2 " + read2
    + " --target " + "flexBarRead1"
    + opts
    + " > " + "flexbarRun1.log")

  if verbose:
    print "running...", flexbar_cmd

  check_call(flexbar_cmd, shell=True)

  #Run flexbar splitslitting on read2
  for r1 in glob.glob("*flexBarRead1*1.fastq"):
    r2 = r1[:-7] + "2.fastq"
    assert os.path.isfile(r1), "read1 of first flexbar run missing"
    assert os.path.isfile(r2), "read2 of first flexbar run missing"

    #Check file has reads in it.
    if os.path.getsize(r1) <= 0:
      print "Skipping empty read files: ", r1[:-7]
      continue


    #Reverse reads so that it splits on the second read file
    flexbar_cmd = (FLEXBAR
      + " --reads " + r2
      + " --reads2 " + r1
      + " --target " + r1[:-7] + "flexBarRead2"
      + opts
      + " > " + r1[:-7] + "_flexbarRun2.log")

    if verbose:
      print "running...", flexbar_cmd

    check_call(flexbar_cmd, shell=True)

  return total_paired_count

def mergePear(pairedMIDs, cpu, verbose):
  #First merge reads
  for r1 in glob.glob("*flexBarRead2*1.fastq"):
    r2 = r1[:-7] + "2.fastq"
    assert os.path.isfile(r1), "read1 of first flexbar run missing"
    assert os.path.isfile(r2), "read2 of first flexbar run missing"

    #Check file has reads in it.
    if os.path.getsize(r1) <= 0:
      print "Skipping empty read files: ", r1[:-7]
      continue

    #The the reads were in the correct orientation orignally so swap back
    if re.match(r'.*forward.*reverse.*', r1):
      r1,r2 = r2,r1

    pear_cmd = (PEAR
      + " -f " + r1
      + " -r " + r2
      + " -o " + r1[:-7] + "PearMerged"
      + " -j " + str(cpu)
      + " -v 10"
      + " -n 100"
      + " > " + r1[:-7] + "pearOut.log")

    if verbose:
      print "Checking read files match..."

    checkFastqPaired(r1, r2, verbose)

    if verbose:
      print "running... ", pear_cmd

    check_call(pear_cmd, shell=True)

  #Now merge files based on MID combinations
  for mid in pairedMIDs:
    filenames = glob.glob("*B_" + mid[0] + "_*B_" + mid[1] + "*.assembled.fastq")
    filenames += glob.glob("*B_" + mid[1] + "_*B_" + mid[0] + "*.assembled.fastq")
    files = []
    for f in filenames:
      #check we have both forward and reverse reads
      if ("_forward" in f) and ("_reverse" in f):
        files.append(f)
    print filenames
    catFiles(files, pairedMIDs[mid] + "_demultiplexTrimMerged.fastq")

  return

def filterReads(filterReads, verbose):
  #Filter errors and convert to fasta
  for f in glob.glob("*_demultiplexTrimMerged.fastq"):

    usearch_filter_cmd = (USEARCH
      + " -fastq_filter"
      + " " + f
      + " -fastaout " + f[:-6] + ".fasta"
      + " -fastqout_discarded " + f[:-6] + "filteredOut.fastq"
      )

    if filterReads:
      usearch_filter_cmd += " -fastq_maxee 1.0"

    usearch_filter_cmd = usearch_filter_cmd + " > " + f[:-6] + "usearchOut.log"

    if verbose:
      print "running... ", usearch_filter_cmd

    check_call(usearch_filter_cmd, shell=True)

def removeLowSupportReads(per_id, min_size, verbose):

  for f in glob.glob("*_demultiplexTrimMerged.fasta"):
    #intially remove singletons
    usearch_cmd = (USEARCH
      + " -derep_prefix"
      + " " + f
      + " -fastaout " + f[:-6] + "_unique.fasta"
      + " -sizeout")

    if verbose:
      print "running... ", usearch_cmd

    check_call(usearch_cmd, shell=True)

    #cluster unique reads and annotate with size
    usearch_cmd = (USEARCH
      + " -cluster_fast"
      + " " + f[:-6] + "_unique.fasta"
      + " -centroids " + f[:-6] + "_centroids.fasta"
      + " -sort size"
      + " -id " + str(per_id))

    if verbose:
      print "running... ", usearch_cmd

    check_call(usearch_cmd, shell=True)

    #map all reads back to otus
    usearch_cmd = (USEARCH
      + " -usearch_global"
      + " " + f
      + " -db " + f[:-6] + "_centroids.fasta"
      + " -strand plus"
      + " -id " + str(per_id)
      + " -dbmatched " + f[:-6] + "_centroids_withSize.fasta"
      + " -sizeout"
      )

    if verbose:
      print "running... ", usearch_cmd

    check_call(usearch_cmd, shell=True)

    #sort by size and filter out low support seqs
    usearch_cmd = (USEARCH
      + " -sortbysize"
      + " " + f[:-6] + "_centroids_withSize.fasta"
      + " -fastaout " + f[:-6] + "_lowSupportFiltered.fasta"
      + " -minsize " + str(min_size)
      )

    if verbose:
      print "running... ", usearch_cmd

    check_call(usearch_cmd, shell=True)

  return

def combineReadFiles(outputfile, verbose):

  if verbose:
    print "combining sample files..."

  #combine read files appending sample name to read headers
  with open(outputfile, 'w') as outfile:
    for f in glob.glob("*_lowSupportFiltered.fasta"):
      for h,s in FastaReader(f):
        outfile.write(">" + h.split("size=")[0]+";sample=" + f[:-59] + "\n")
        outfile.write(s + "\n")

  return

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
  with open("hmmerDBLalphaSearch.txt", 'r') as infile:
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

def countReads(fastafile):
  count = 0
  for h,s in FastaReader(fastafile):
    count += 1
  return count

def readPearLog(filename):
  pearResults = {}
  with open(filename, 'r') as infile:
    for line in infile:
      if "Assembled reads ..." in line:
        pearResults["total"] = int(line.split()[5].replace(",",""))
        pearResults["assembled"] = int(line.split()[3].replace(",",""))
      if "Discarded reads .." in line:
        pearResults["discarded"] = int(line.split()[3].replace(",",""))
      if "Not assembled reads" in line:
        pearResults["not_assembled"] = int(line.split()[4].replace(",",""))

  if len(pearResults.keys())!=4:
    MyError("Problem reading in Pear log file: " + filename)

  return pearResults

def readFlexBarLog(filename):
  flexbarResults = {}
  with open(filename, 'r') as infile:
    AtStats = False
    for line in infile:
      line = line.strip()
      if AtStats:
        if "Processed reads" in line:
          flexbarResults["processed_reads"] = int(line.split()[2])
        elif "skipped due to uncalled bases" in line:
          flexbarResults["skip_uncalled_bases"] = int(line.split()[5])
        elif "skipped unassigned reads" in line:
          flexbarResults["skip_unassigned_reads"] = int(line.split()[3])
        elif "finally skipped short reads" in line:
          flexbarResults["skip_short_reads"] = int(line.split()[4])
        elif "skipped single paired reads" in line:
          flexbarResults["skip_single_reads"] = int(line.split()[4])
        elif "Discarded reads overall" in line:
          flexbarResults["discarded_reads"] = int(line.split()[3])
      elif "Filtering statistics" in line:
        AtStats=True
        infile.next()

  if len(flexbarResults.keys())!=6:
    MyError("Problem reading flexbar log: " + filename)

  return flexbarResults

def calculateSummaryStatistics(pairedMIDs, outputfile, total_reads
  , total_reads_before_contaminant_filtering
  , total_reads_after_contaminant_filtering
  , verbose):

  ##Calculate flexbar summary statisics
  flexBarStats = readFlexBarLog("flexbarRun1.log")
  logfiles = glob.glob("*flexbarRun2.log")
  for f in logfiles:
    tempLog = readFlexBarLog(f)
    flexBarStats["processed_reads"] += tempLog["processed_reads"]
    flexBarStats["skip_uncalled_bases"] += tempLog["skip_uncalled_bases"]
    flexBarStats["skip_unassigned_reads"] += tempLog["skip_unassigned_reads"]
    flexBarStats["skip_short_reads"] += tempLog["skip_short_reads"]
    flexBarStats["skip_single_reads"] += tempLog["skip_single_reads"]
    flexBarStats["discarded_reads"] += tempLog["discarded_reads"]


  ##Calculate pear merging statistics
  #Now merge files based on MID combinations
  pearStats = {}
  for mid in pairedMIDs:
    filenames = glob.glob("*B_" + mid[0] + "_*B_" + mid[1] + "*pearOut.log")
    filenames += glob.glob("*B_" + mid[1] + "_*B_" + mid[0] + "*pearOut.log")
    sample = pairedMIDs[mid]
    pearStats[sample] = readPearLog(filenames[0])
    for f in filenames[1:]:
      tempPear = readPearLog(f)
      pearStats[sample]["total"] += tempPear["total"]
      pearStats[sample]["assembled"] += tempPear["assembled"]
      pearStats[sample]["discarded"] += tempPear["discarded"]
      pearStats[sample]["not_assembled"] += tempPear["not_assembled"]

  ##Calculate support filtering statistics
  #Count number of reads before support filtering
  initalCounts = {}
  for f in glob.glob("*_demultiplexTrimMerged.fasta"):
    sample = f.split("_demultiplexTrimMerged")[0]
    initalCounts[sample] = countReads(f)

  #Count number of centroids
  centroidCounts = {}
  for f in glob.glob("*_demultiplexTrimMerged_centroids.fasta"):
    sample = f.split("_demultiplexTrimMerged")[0]
    centroidCounts[sample] = countReads(f)

  #Count number of centroids with sufficient support
  supportCentroids = {}
  for f in glob.glob("*_demultiplexTrimMerged_lowSupportFiltered.fasta"):
    sample = f.split("_demultiplexTrimMerged")[0]
    supportCentroids[sample] = countReads(f)

  ##write to outputfile
  with open(outputfile, 'w') as outfile:
    outfile.write("# Total paired end reads: " + str(total_reads) + "\n")
    outfile.write("# Total reads prior to contaminant filtering: "
      + str(total_reads_before_contaminant_filtering) + "\n")
    outfile.write("# Total reads after contaminant filtering: "
      + str(total_reads_after_contaminant_filtering) + "\n")

    outfile.write("#######  FlexBar Summary statistics  #######\n")
    outfile.write("# Total reads processed: " + str(flexBarStats["processed_reads"]) + "\n")
    outfile.write("# Reads skipped due to uncalled bases: " + str(flexBarStats["skip_uncalled_bases"]) + "\n")
    outfile.write("# Unassigned reads: " + str(flexBarStats["skip_unassigned_reads"]) + "\n")
    outfile.write("# Short reads skipped: " + str(flexBarStats["skip_short_reads"]) + "\n")
    outfile.write("# Single reads skipped: " + str(flexBarStats["skip_single_reads"]) + "\n")
    outfile.write("# Total reads discarded: " + str(flexBarStats["discarded_reads"]) + "\n")
    p = 100.0*(1-flexBarStats["discarded_reads"]/float(flexBarStats["processed_reads"]))
    outfile.write("# Proportion of reads kept: " + str(p) + "\n")

    outfile.write("#######  Sample specific summary statistics  #######\n")
    #now write out sample specific statistics
    outfile.write(",".join(["Sample", "PreMerge", "Merged", "Discarded"
      , "Not Assembled", "Filtered", "Centroids", "Centroids with support"])+"\n")

    for sample in pearStats:
      outfile.write(",".join([sample
        , str(pearStats[sample]["total"])
        , str(pearStats[sample]["assembled"])
        , str(pearStats[sample]["discarded"])
        , str(pearStats[sample]["not_assembled"])
        , str(pearStats[sample]["total"] - initalCounts[sample])
        , str(centroidCounts[sample])
        , str(supportCentroids[sample])
        ])+ "\n")

  return


def main():
  parser = argparse.ArgumentParser(description='Process paired end Illumina reads to obtain DBLalpha OTUs.')

  parser.add_argument('-o', '--outputDir', action=writeable_dir
    , dest='outputdir'
    , help="location of output directory. Will be created if it doesn't exist"
    , required=True)

  parser.add_argument('-r', '--read1', dest='read1', action=readable_file
    , help="location of first read file."
    , required=True)

  parser.add_argument('-R', '--read2', dest='read2', action=readable_file
    , help="location of paired read file."
    , required=True)

  parser.add_argument('-d', '--desc', dest='desc', action=readable_file
    , help="location of desc mapping file."
    , required=True)

  parser.add_argument('--filter', dest='filterReads', action='store_true'
    , default=False
    , help='filter merged reads using usearch if more than 1 expected error (default=False)')

  parser.add_argument('--minSize', dest='min_size', type=int, default=4
    , help="minimum support for a read to be kept. (default=4)")

  parser.add_argument('--perID', dest='perID', type=float, default=0.96
    , help="percentage ID threshold. (default=0.96)")

  parser.add_argument('--cpu', dest='cpu', type=int, default=1
    , help="number of cpus to use. (default=1)")

  parser.add_argument('--verbose', dest='verbose', action='store_true'
    , default=False
    , help='print verbose output (default=False)')

  args = parser.parse_args()

  #get full path names
  args.outputdir = os.path.abspath(args.outputdir) + "/"
  args.read1 = os.path.abspath(args.read1)
  args.read2 = os.path.abspath(args.read2)
  args.desc = os.path.abspath(args.desc)

  #create the output directory if it doesn't exist
  try:
      os.mkdir(args.outputdir)
  except OSError, e:
      if e.errno != 17: #ignores error if folders has already been created.
          raise
      pass

  #first create a temporary directory in the output directory
  try:
      outputdir = args.outputdir + "temp_files/"
      os.mkdir(outputdir)
  except OSError, e:
      if e.errno != 17: #ignores error if folders has already been created.
          raise
      pass

  curr_dir = os.getcwd()
  os.chdir(outputdir)

  barcodefile, pairedMIDs = convertDescToFlexBarcodes(args.desc
    , outputdir + "barcodesWithPrimer.fasta"
    , args.verbose)

  total_paired_count = demultiplexTrim(barcodefile, args.read1, args.read2
    , args.cpu, args.verbose)

  mergePear(pairedMIDs, args.cpu, args.verbose)

  filterReads(args.filterReads, args.verbose)

  removeLowSupportReads(args.perID, args.min_size, args.verbose)

  combinedFile = os.path.basename(args.read1[:-7])+"_combinedNotClean.fasta"
  combineReadFiles(combinedFile, args.verbose)

  filterWithHMMER(combinedFile
    , args.outputdir + os.path.basename(args.read1[:-7])
    , args.cpu, args.verbose)

  #Get contaminant filtering counts
  total_reads_before_contaminant_filtering = countReads(combinedFile)
  total_reads_after_contaminant_filtering = countReads(args.outputdir
    + os.path.basename(args.read1[:-7]) + "_DBLa_cleaned.fasta")

  calculateSummaryStatistics(pairedMIDs
    , args.outputdir + "summaryStatistics.log"
    , total_paired_count
    , total_reads_before_contaminant_filtering
    , total_reads_after_contaminant_filtering
    , args.verbose)



if __name__ == '__main__':
  main()








