#!/usr/bin/python
import os
import sys
import glob
import string
import operator 
from multiprocessing import Pool, Process
from string import atof
from itertools import imap
from Bio import SeqIO

def rc(seq):
  complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
  rcseq = seq.translate(complements)[::-1]
  return rcseq

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(imap(operator.ne, str1, str2))

def offsetcheck(R_seq,refseqs):
  for i in range(0,len(R_seq)-240+1):
    Fseq = R_seq[i:i+240]
    Rseq = rc(Fseq)
    for ref in refseqs.keys():
      refseq  = refseqs[ref]
      Fhdist  = hamming(refseq,Fseq)
      Rhdist  = hamming(refseq,Rseq)
      if Fhdist <= 6: 
        strand = 'F'
        Roffset = i
        return [ref,Roffset,strand]
      if Rhdist <= 6: 
        strand = 'R'
        Roffset = i
        return [ref,Roffset,strand]
  return 'bad'

def MapNPair(R1file, R2file,mfile,refseqs):
  workID = R1file
  print 'working on %s' % workID
  R1file = SeqIO.parse(R1file, "fastq")
  R2file = SeqIO.parse(R2file, "fastq")
  mfile  = open(mfile,'w')
  badBC = 0
  total = 0
  for R1record in R1file:
    R2record = R2file.next()
    R1_ID    = str(R1record.id)
    R1_bc    = str(R1record.seq)[0:4]
    R1_seq   = str(R1record.seq)[4::]
    R1_qual  = map(int,R1record.letter_annotations["phred_quality"])
    R2_ID    = str(R2record.id)
    R2_bc    = str(R2record.seq)[0:4]
    R2_seq   = str(R2record.seq)[4::]
    R2_qual  = map(int,R2record.letter_annotations["phred_quality"])
    #QUALITY CONTROL#
    assert(R1_ID == R2_ID)
    total += 1
    if R1_bc != R2_bc: 
      badBC += 1
      continue
    if len(R1_seq) < 240: continue
    if len(R2_seq) < 240: continue
    #END OF QUALITY CONTROL#
    #EXTRACT OFFSET INFO
    R1_info   = offsetcheck(R1_seq,refseqs)
    R2_info   = offsetcheck(R2_seq,refseqs)
    #QUALITY CONTROL#
    if R1_info == 'bad' or R2_info == 'bad': continue
    if R1_info[0] != R2_info[0]: continue
    if R1_info[2] == R2_info[2]: continue
    #END OF QUALITY CONTROL#
    #CALL MUTATION
    WT_Amp    = R1_info[0]
    refseq    = refseqs[WT_Amp]
    R1_offset = R1_info[1]
    R1_strand = R1_info[2]
    R2_offset = R2_info[1]
    R2_strand = R2_info[2]
    R1_seq    = R1_seq[R1_offset:R1_offset+240]
    R1_qual   = R1_qual[R1_offset:R1_offset+240]
    R2_seq    = R2_seq[R2_offset:R2_offset+240]
    R2_qual   = R2_qual[R2_offset:R2_offset+240]
    if R1_strand == 'R': R1_seq = rc(R1_seq); R1_qual.reverse()
    if R2_strand == 'R': R2_seq = rc(R2_seq); R2_qual.reverse()
    if R1_seq != R2_seq: continue
    Muts = []
    for n in range(0,len(refseq)):
      if R1_seq[n] != refseq[n] and R1_seq[n] == R2_seq[n]:
        Mut = refseq[n]+str(n+1)+R1_seq[n]
        Muts.append(Mut)
    if len(Muts) == 0:
      Muts = ['WT']
    mfile.write(WT_Amp+"\t"+R1_bc[0:3]+"\t"+'-'.join(Muts)+"\n")
  R1file.close()
  R2file.close()
  mfile.close()
  print 'For', workID, badBC, 'bad barcodes out of:', total


#READ IN REFERENCE SEQUENCE

def hashinref(reffile):
  reffile = open(reffile,'r')
  refseqs = {}
  for line in reffile.xreadlines():
    if '>' in line: 
      ID = line.rstrip().replace('>','')
    else:
      refseqs[ID] = line.rstrip()
  return refseqs


#WRAPPER
def multi_run_wrapper(args):
   return MapNPair(*args)


#MAIN#
def main():
  refseqs = hashinref('Fasta/flu1amp.fa')
  R1files = sorted(glob.glob('fastq/*_R1_*.fastq'))
  R2files = [R1file.replace('_R1_','_R2_') for R1file in R1files]
  fileIDs = [R1file.rsplit('_')[4].rsplit('.')[0] for R1file in R1files]
  mfiles  = ['tmp/'+fileID+'.m' for fileID in fileIDs]
  pool   = Pool(processes=2)
  Inputs = []
  [Inputs.append((R1files[n],R2files[n],mfiles[n],refseqs)) for n in range(len(R1files))]
  pool.map(multi_run_wrapper,Inputs)
  pool.close()
  pool.join()
  os.system('cat tmp/*.m > result/AllM_1')
  os.system('rm tmp/*')
  #os.system('play --no-show-progress --null --channels 1 synth %s sine %f' % ( 0.5, 2000))

if __name__ == "__main__":
  main()
