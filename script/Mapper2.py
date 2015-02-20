#!/usr/bin/python
import os
import sys
import glob

def readframe(aa, pos,firstbaseORF):
  wtaa  = aa[0]
  mutaa = aa[1]
  pos   = (pos-firstbaseORF)/3+1
  return wtaa+str(pos)+mutaa

def adjustmutpos(muts, offset):
  if muts == 'WT': return muts
  else:
    muts  = muts.rsplit('-')
    mlist = []
    for m in muts:
      pos = str(int(m[1:-1])+offset-1)
      mlist.append(m[0]+pos+m[-1])
    return '-'.join(mlist)

#READ IN OFFSET FILE
def hashinoffset(filename):
  infile  = open(filename,'r')
  offsetH = {}
  for line in infile.xreadlines():
    line = line.rstrip().rsplit("\t")
    offsetH[line[0]] = [int(line[1]),int(line[2])]
  infile.close()
  return offsetH

#READ IN BARCODE FILE
def hashinBC(filename):
  infile   = open(filename,'r')
  barcodes = {}
  pops     = []
  for line in infile.xreadlines():
    line = line.rstrip().rsplit("\t")
    barcodes[line[0]] = line[1]
    pops.append(line[1])
  infile.close()
  return pops, barcodes

#READ IN PROTEIN LEVEL DATA
def hashinaa(filename,seg):
  infile = open('Fasta/fluTranslate','r')
  Rframe = {}
  for line in infile.xreadlines():
    array = line.rstrip().rsplit("\t")
    ID = ''.join(array[0:3])
    ref = array[3]
    WTaa1 = array[4]
    Mutaa1 = array[5]
    WTaa2 = array[6]
    Mutaa2 = array[7]
    if ref == seg: Rframe[ID] = [WTaa1+Mutaa1, WTaa2+Mutaa2]
  infile.close()
  return Rframe

#GENOTYPE AND DEPTH HASH INITIATE
def varinitiate(pops,offsetH):
  GenotypeH = {}
  depthH    = {}
  WTcount   = {}
  for pop in pops:
    GenotypeH[pop] = {}
    depthH[pop]    = {}
    WTcount[pop]   = {}
    for amp in offsetH.keys():
      depthH[pop][amp] = 0
      WTcount[pop][amp] = 0
  return GenotypeH, depthH, WTcount

#Hashing in Genotypes
def hashinG(filename,GenotypeH,depthH,WTcount,barcodes,offsetH):
  infile = open(filename,'r')
  for line in infile.xreadlines():
    line = line.rstrip().rsplit("\t")
    if line[1] not in barcodes.keys(): continue
    amp    = line[0]
    offset = offsetH[amp][0]
    pop    = barcodes[line[1]]
    muts   = line[2]
    muts   = adjustmutpos(muts,offset)
    if 'C682G' in muts: print muts, t
    depthH[pop][amp] += 1
    if muts == 'WT':
      WTcount[pop][amp] += 1
    else:
      if GenotypeH[pop].has_key(muts): 
        GenotypeH[pop][muts] += 1
      else: 
        GenotypeH[pop][muts] = 1
  infile.close()
  return GenotypeH,depthH

#SUMMARIZE GENOTYPE
def UniqGenotype(GenotypeH,pops):
  Genotypes = []
  for pop in pops:
    Genotypes.extend(GenotypeH[pop].keys())
  Genotypes = list(set(Genotypes))
  return Genotypes

#WTCOUNT AND DEPTH OUTPUT
def WTnDepthout(filename,WTcount,depthH,pops,offsetH):
  outfile = open(filename,'w')
  outfile.write('pop'+"\t"+'amp'+"\t"+'WTcount'+"\t"+'Depth'+"\n")
  for pop in pops:
    for amp in offsetH.keys():
      outfile.write("\t".join([pop, amp, str(WTcount[pop][amp]), str(depthH[pop][amp])])+"\n")
  outfile.close()

#GENOTYPES
def GenoOut(outfileA,outfileS,pops,GenotypeH,Genotypes):
  outfileA = open(outfileA,'w')
  outfileS = open(outfileS,'w')
  header   = 'Genotype'+"\t"+"\t".join(pops)
  outfileA.write(header+"\n")
  outfileS.write(header+"\n")
  for G in Genotypes:
    counts = []
    for pop in pops:
      if GenotypeH[pop].has_key(G): 
        counts.append(GenotypeH[pop][G])
      else: 
        counts.append(0)
    out = G+"\t"+"\t".join(map(str,counts))
    outfileA.write(out+"\n")
    if '-' not in G:
      outfileS.write(out+"\n")
  outfileA.close()
  outfileS.close()

#Single mutation count writing file
def SmutOut(filename,pops,Genotypes,Rframe,offsetH,WTcount,depthH,GenotypeH,firstbaseORF):
  outfileC = open(filename,'w')
  header   = ['Pos','Genotype','Frame1','Frame2']
  for pop in pops:
    header.extend([pop, pop+'_WT',pop+'_Dep'])
  outfileC.write("\t".join(header)+"\n")

  for G in Genotypes:
    if '-' not in G and 'N' not in G:
      pos = int(G[1:-1])
      F1  = readframe(Rframe[G][0],pos,firstbaseORF)
      F2  = readframe(Rframe[G][1],8,firstbaseORF)
      out = [str(pos),G,F1,F2]
      for pop in pops:
        wtc = 0
        dep = 0 
        for amp in offsetH.keys():
          if pos >= offsetH[amp][0] and pos <= offsetH[amp][1]:
            wtc += WTcount[pop][amp]
            dep += depthH[pop][amp]
        if GenotypeH[pop].has_key(G): out.append(str(GenotypeH[pop][G]))
        else: out.append(str(0))
        out.append(str(wtc))
        out.append(str(dep))
      outfileC.write("\t".join(out)+"\n")
  outfileC.close()
  
def main():
  firstbaseORF   = 16 #first base for the ORF (position of the first nucleotide coding for the starting methionine)
  seg            = 'flu1'
  fluinfofile    = 'Fasta/fluTranslate'
  AllM           = 'result/AllM_1' #Output file from Mapper1
  offsetfile     = 'Fasta/flu1offset'
  barcodefile    = 'Fasta/BarCode'
  WTnDepthOfile  = 'result/WT_N_Depth'
  AllGenoOfile   = 'result/AGenotypes'
  MonoGenoOfile  = 'result/SGenotypes'
  SmutOfile      = 'result/SMut'
  Rframe         = hashinaa(fluinfofile,seg)
  offsetH        = hashinoffset(offsetfile)
  pops, barcodes = hashinBC(barcodefile)
  GenotypeH, depthH, WTcount = varinitiate(pops,offsetH)
  GenotypeH, depthH          = hashinG(AllM,GenotypeH,depthH,WTcount,barcodes,offsetH)
  Genotypes                  = UniqGenotype(GenotypeH,pops)
  WTnDepthout(WTnDepthOfile,WTcount,depthH,pops,offsetH)
  GenoOut(AllGenoOfile,MonoGenoOfile,pops,GenotypeH,Genotypes)
  SmutOut(SmutOfile,pops,Genotypes,Rframe,offsetH,WTcount,depthH,GenotypeH,firstbaseORF)

if __name__ == "__main__":
  main()

