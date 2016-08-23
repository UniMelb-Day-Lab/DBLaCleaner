from mungo.fasta import FastaReader
import os,sys
import glob


for iso in ["3D7", "[P1]HB3.MID", "HB3xDD2", "DD2xHB3", "[P1]DD2.MID"]:
  with open("combined_"+iso+".fasta", 'w') as outfile:
    for f in glob.glob("./chimeric_filtered*/*"+iso+"*ncdenovo.fasta"):
        for h,s in FastaReader(f):
          sample = f.split("/")[-1].split("_demultiplexTrimMerged_")[0]
          outfile.write(">" + h + ";sample=" + sample + "_ncdenovo;\n")
          outfile.write(s + "\n")
    for f in glob.glob("./chimeric_filtered*/*"+iso+"*_otus.fa"):
        for h,s in FastaReader(f):
          sample = f.split("/")[-1].split("_demultiplexTrimMerged_")[0]
          outfile.write(">" + h + ";sample=" + sample + "_otuclust;\n")
          outfile.write(s + "\n")
    for f in glob.glob("./chimeric_filtered*/*"+iso+"*ncref.fasta"):
        for h,s in FastaReader(f):
          sample = f.split("/")[-1].split("_demultiplexTrimMerged_")[0]
          outfile.write(">" + h + ";sample=" + sample + "_ncref;\n")
          outfile.write(s + "\n")
    for f in glob.glob("./*"+iso+"*DBLa_cleaned.fasta"):
        for h,s in FastaReader(f):
          sample = f.split("/")[-1].split("_demultiplexTrimMerged_")[0]
          outfile.write(">" + h + ";sample=" + sample + ";\n")
          outfile.write(s + "\n")

