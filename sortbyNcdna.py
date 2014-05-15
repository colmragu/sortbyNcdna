from __future__ import print_function
from Bio import SeqIO
import pdb
import sys

def sortbyNcdna(sortbyNCDNA, gb_file):
  scan_pos = 0
  f = open(gb_file,"r")
  ncdna_last = "S"

  count =0
  try:
    gb_record = SeqIO.read(f, "genbank")
  except:
    print ("Skipping %s , cannot read file"%gb_file)
    return (sortbyNCDNA) 
  f.close()
  for feature in gb_record.features[1:]: #Select all but the first entry in gbk file
    feature = pad_gb_features(feature)
    ncdna = str(gb_record.seq[scan_pos:feature.location.start.position+1]) #ncdna = non coding dna
    if ncdna == "":
      sortbyNCDNA = addfeaturetodict(ncdna_last,sortbyNCDNA, feature, gb_file)
    else:
      sortbyNCDNA = addfeaturetodict(ncdna,sortbyNCDNA, feature, gb_file)
      ncdna_last= ncdna
    scan_pos = feature.location.end.position
  return(sortbyNCDNA)

def addfeaturetodict(ncdna,sortbyNCDNA, feature, organism):
  if not sortbyNCDNA.has_key(ncdna):
    sortbyNCDNA[ncdna]={}
  if not sortbyNCDNA[ncdna].has_key(feature.qualifiers["db_xref"][0]):
    sortbyNCDNA[ncdna][feature.qualifiers["db_xref"][0]] ={} 

  #sortbyNCDNA[ncdna][feature.qualifiers["db_xref"][0]] = feature.qualifiers["note"].join(feature.location.start).join(feature.location.end)
  sortbyNCDNA[ncdna][feature.qualifiers["db_xref"][0]] = feature.qualifiers["note"]  ,int (feature.location.start) ,int (feature.location.end) ,organism #add more information

  return(sortbyNCDNA) 

def pad_gb_features(feature):
  ### Make sure all gb objects conatain a dic entry for note, db_xref
  if not feature.qualifiers.has_key("note"):
    feature.qualifiers["note"] = [""]
  if not feature.qualifiers.has_key("db_xref"):
    feature.qualifiers["db_xref"]=["nodb_xref"]
  return feature

def print_sorted(ncdna_sort, outputfile):
  f = open(outputfile, "w")
  for ncdna in sorted(ncdna_sort.keys(), key=len): 
    print (ncdna, file=f)
    for gene_id in ncdna_sort[ncdna]:
      print (gene_id, file=f, end=" ")

      print (ncdna_sort[ncdna][gene_id], file=f)
    print ("", file=f)
  f.close()

def sortall(gb_files): 
  ncdna_sort = {}
  for gb_file in gb_files:
    print(gb_file)
    ncdna_sort = sortbyNcdna(ncdna_sort, gb_file)

  return(ncdna_sort)

def parse_commandline():
  help_message = "Incorrect command you need to supply a gbk file \n python sortbyNcdna.py data/NC_009930.gbk outputfile \n python sortbyNcdna.py data/NC_017100.gbk data/NC_017118.gbk data/NC_009931.gbk outputfile \n python sortbyNcdna.py data/*.gbk outputfile"
  if len(sys.argv) > 1:
    gb_files = sys.argv[1:-1]  
    outputfile = sys.argv[-1]
  else:
    print(help_message)
    return(0)
  ncdna_sorted = sortall(gb_files)

  print_sorted(ncdna_sorted, outputfile)

parse_commandline()
