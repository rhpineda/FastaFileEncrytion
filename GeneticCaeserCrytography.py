#-----------------------
#NOTES
# 1 is file
# Caeser cipher fxn
# Transcription fxn
# -include B/Z?
# -include O?
# -punctuation? (-)
# Translation fxn
# -different genetic codes, default is canonical
# reverse fxns too
#A[B]CDEFGHI((J))KLMN[O]PQRST((U))VW((X))Y[Z]
#B (AAT) from N 
#Z (CAA) from Q
#O {3 codons} -1 for period
#period?? (get amber stop?)
#J AAA
#U CCC
#X GGG
#space TTT
#??? GTG

#maybe do lists of lists to test speed

#test -[caeser]-> uftu -[tl]-> TACGTGAC -[tl]-> uftu -[caeser]-> test 

'''
GOAL
1. Make a fake fasta file that encodes some text
  -Encode in some fake information
  -have a way to allow for multiple sequences for multiple files?
  -randomized caeser +way to decode
2. Have realistic aa freq
  -randomize which codon to use for ones with multiple codons
3. have some way to encode some random seq "gene" that says how to decode
4. expand more ascii
5. fna, faa, frn, fa
'''
#-----------------------------
#CODE
import sys
import gzip
import random
import argparse

parser = argparse.ArgumentParser(description= "Program to encode txt -> .fna")
parser.add_argument('-s', required=True, type=str, 
	metavar='<str>', help='file to evaluate')
parser.add_argument('-c', required=False, type=int, default=0,
	metavar='<int>', help='"Caeser Shift Amount", default is: 0')
parser.add_argument('--decode', action = 'store_true', \
  help='"Encode/Decode", default is: True, Encode')
arg = parser.parse_args()

def caeser(seq, shift = 0):
  seq = seq.upper()
  caeserdict1 = { 
  'A':1, 'B':2, 'C':3, 'D':4, 'E':5, 'F':6, 'G':7, 'H':8, 'I':9, 'J':10,
  'K':11, 'L':12, 'M':13, 'N':14, 'O':15, 'P':16, 'Q':17, 'R':18, 'S':19,
  'T':20, 'U':21, 'V':22, 'W':23, 'X':24, 'Y':25, 'Z':26, '.':27, ' ':28,
  '#':29
  } #no shift
  caeserdict2 = { 
  'A':1, 'B':2, 'C':3, 'D':4, 'E':5, 'F':6, 'G':7, 'H':8, 'I':9, 'J':10,
  'K':11, 'L':12, 'M':13, 'N':14, 'O':15, 'P':16, 'Q':17, 'R':18, 'S':19,
  'T':20, 'U':21, 'V':22, 'W':23, 'X':24, 'Y':25, 'Z':26, '.':27, ' ':28,
  '#':29
  } #shifts
  newdict = {} #used in encryption
  encodedseq = ''  #working new seq
  for i in caeserdict1: 
    caeserdict2[i] = (caeserdict1[i]+shift)%29 # shift caeser dict 2
  for key1, value1 in caeserdict1.items():
    for key2, value2 in caeserdict2.items():
      if value1 == value2:
        newdict.update({key1:key2}) #makes newdict
  for i in range(len(seq)):
    if seq[i] in newdict:
      encodedseq += newdict[seq[i]]
    else:
      encodedseq += '#'
  return(encodedseq)
def tl(seq, encode = True): #encode/translate
  endict  = {
	'A':'GCT', 'A':'GCC', 'A':'GCA', 'A':'GCG',
	'B':'AAT', 'C':'TGT', 'C':'TGC', 'D':'GAT',
	'D':'GAC', 'E':'GAA', 'E':'GAG', 'F':'TTC',
	'G':'GGT', 'G':'GGC', 'G':'GGA', 'H':'CAT', 
	'H':'CAC', 'I':'ATT', 'I':'ATC', 'I':'ATA', 
	'J':'AAA', 'K':'AAG', 'L':'TTA', 'L':'TTG',
	'L':'CTT', 'L':'CTC', 'L':'CTA', 'L':'CTG',
	'M':'ATG', 'N':'AAC', 'O':'TAA', 'O':'TGA', 
	'P':'CCT', 'P':'CCA', 'P':'CCG', 'Q':'CAG',
  'R':'CGT', 'R':'CGC', 'R':'CGA', 'R':'CGG', 
  'R':'AGA', 'R':'AGG', 'S':'TCT', 'S':'TCC',
  'S':'TCA', 'S':'TCG', 'S':'AGT', 'S':'AGC',
  'T':'ACT', 'T':'ACC', 'T':'ACA', 'T':'ACG',
  'U':'CCC', 'V':'GTT', 'V':'GTC', 'V':'GTA',
  'W':'TGG', 'X':'GGG', 'Y':'TAT', 'Y':'TAC', 
  'Z':'CAA', '.':'TAG', ' ':'TTT', '#':'GTG'
  }
  decdict = {
	'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
	'AAT':'B', 'TGT':'C', 'TGC':'C', 'GAT':'D',
	'GAC':'D', 'GAA':'E', 'GAG':'E', 'TTC':'F',
	'GGT':'G', 'GGC':'G', 'GGA':'G', 'CAT':'H', 
	'CAC':'H', 'ATT':'I', 'ATC':'I', 'ATA':'I', 
	'AAA':'J', 'AAG':'K', 'TTA':'L', 'TTG':'L',
	'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
	'ATG':'M', 'AAC':'N', 'TAA':'O', 'TGA':'O', 
	'CCT':'P', 'CCA':'P', 'CCG':'P', 'CAG':'Q',
  'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 
  'AGA':'R', 'AGG':'R', 'TCT':'S', 'TCC':'S',
  'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
  'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
  'CCC':'U', 'GTT':'V', 'GTC':'V', 'GTA':'V',
  'TGG':'W', 'GGG':'X', 'TAT':'Y', 'TAC':'Y', 
  'CAA':'Z', 'TAG':'.', 'TTT':' ', 'GTG':'#'
  }
  seq = seq.upper()
  encodedseq = ''
  
  if encode == True:
    for i in range(len(seq)):
      if seq[i] in endict:
        encodedseq += endict[seq[i]]
      else:
        encodedseq += 'GTG'
    return(encodedseq)    
  else:
    for i in range(0, len(seq), 3):
      if str(seq[i:i+3]) in decdict:
        encodedseq += decdict[str(seq[i:i+3])]
    return(encodedseq)
def accession():
  accession = ''.join(random.choice('ABCDEFGHIJKLMNOPQRSTUVWZYZ') \
  for i in range(random.randint(1,3)))
  accession += str(random.randint(10000,99999))
  if (random.random()> 0.7): accession += ('.' + str(random.randint(1,9)))
  return(accession)
def header():
  dba = ['gb', 'emb', 'pir', 'sp', 'ref', 'dbj', 'prf', 'tpg', 'tpe', 'tpd', 'tr']
  info1 = ['Homo sapiens', 'A.thaliana', 'D.melanogaster', 'M.musculus', 
            'Xenopus laevis', 'C.elegans', 'Saccharomyces cerevisiae']
  info2 = ['putative transcription factor', 'norovirus', 'incomplete genome', 
            'complete genome', "RNA polymerase III 5'UTR", 
            "RNA polymerase 5'UTR", "unknown nuclease 3'UTR",
            'partial genome']
  #Generate Accession
  
  #generate Info Line
  info  = (
    '>' + random.choice(dba) + '|'  + accession() + '|' +
    random.choice(info1) + '|' +  random.choice(info2) + '\n'
    )
  return(info)

def mff(seq, shift = 0):
  enc = ''
  encw = enc
  enc += tl(caeser(seq,0), True)
  for i in range(0, len(enc), 60):
      encw += enc[i:i+60] + '\n'
  title = accession() + '00000' + str(random.randint(10000,99999)) + '.fna'
  with open(title,'w') as file:
    file.write(header() + encw)
def rff(seq, shift = 0):
  enc = ''
  encr = enc
  encr += tl(seq, False)
  with open('decoded.txt','w') as file:
    file.write(encr)

if arg.decode == False:
  with open(arg.s, 'r') as file:
    file_contents = file.read()
  mff(file_contents)
  print("File Created")
else:
  with open(arg.s, 'rt') as fp:
    seq = ''
    for line in fp.readlines():
      if not line.startswith('>'):
        seq += line
        seq = seq.replace('\n','')
  rff(seq)
  print("File Decoded")
#-------------------------------------------------
#unused code
'''
def enctl(seq): #encode/translate
  translatedict  = {
	'A':'GCT', 'A':'GCC', 'A':'GCA', 'A':'GCG',
	'B':'AAT', 'C':'TGT', 'C':'TGC', 'D':'GAT',
	'D':'GAC', 'E':'GAA', 'E':'GAG', 'F':'TTC',
	'G':'GGT', 'G':'GGC', 'G':'GGA', 'H':'CAT', 
	'H':'CAC', 'I':'ATT', 'I':'ATC', 'I':'ATA', 
	'J':'AAA', 'K':'AAG', 'L':'TTA', 'L':'TTG',
	'L':'CTT', 'L':'CTC', 'L':'CTA', 'L':'CTG',
	'M':'ATG', 'N':'AAC', 'O':'TAA', 'O':'TGA', 
	'P':'CCT', 'P':'CCA', 'P':'CCG', 'Q':'CAG',
  'R':'CGT', 'R':'CGC', 'R':'CGA', 'R':'CGG', 
  'R':'AGA', 'R':'AGG', 'S':'TCT', 'S':'TCC',
  'S':'TCA', 'S':'TCG', 'S':'AGT', 'S':'AGC',
  'T':'ACT', 'T':'ACC', 'T':'ACA', 'T':'ACG',
  'U':'CCC', 'V':'GTT', 'V':'GTC', 'V':'GTA',
  'W':'TGG', 'X':'GGG', 'Y':'TAT', 'Y':'TAC', 
  'Z':'CAA', '.':'TAG', ' ':'TTT', '#':'GTG'
  }
  seq = seq.upper()
  encodedseq = ''
  
  for i in range(len(seq)):
    if seq[i] in translatedict:
      encodedseq += translatedict[seq[i]]
    else:
      encodedseq += 'GTG'
  return(encodedseq)
def dectl(seq):
  translatedict  = {
	'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
	'AAT':'B', 'TGT':'C', 'TGC':'C', 'GAT':'D',
	'GAC':'D', 'GAA':'E', 'GAG':'E', 'TTC':'F',
	'GGT':'G', 'GGC':'G', 'GGA':'G', 'CAT':'H', 
	'CAC':'H', 'ATT':'I', 'ATC':'I', 'ATA':'I', 
	'AAA':'J', 'AAG':'K', 'TTA':'L', 'TTG':'L',
	'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
	'ATG':'M', 'AAC':'N', 'TAA':'O', 'TGA':'O', 
	'CCT':'P', 'CCA':'P', 'CCG':'P', 'CAG':'Q',
  'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 
  'AGA':'R', 'AGG':'R', 'TCT':'S', 'TCC':'S',
  'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
  'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
  'CCC':'U', 'GTT':'V', 'GTC':'V', 'GTA':'V',
  'TGG':'W', 'GGG':'X', 'TAT':'Y', 'TAC':'Y', 
  'CAA':'Z', 'TAG':'.', 'TTT':' ', 'GTG':'#'
  }
  seq = seq.upper()
  encodedseq = ''
  for i in range(0, len(seq), 3):
    if str(seq[i:i+3]) in translatedict:
      encodedseq += translatedict[str(seq[i:i+3])]
  return(encodedseq)
translatedict  = {
	'A':'GCT', 'A':'GCC', 'A':'GCA', 'A':'GCG',
	'B':'AAT', 'C':'TGT', 'C':'TGC', 'D':'GAT',
	'D':'GAC', 'E':'GAA', 'E':'GAG', 'F':'TTC',
	'G':'GGT', 'G':'GGC', 'G':'GGA', 'H':'CAT', 
	'H':'CAC', 'I':'ATT', 'I':'ATC', 'I':'ATA', 
	'J':'AAA', 'K':'AAG', 'L':'TTA', 'L':'TTG',
	'L':'CTT', 'L':'CTC', 'L':'CTA', 'L':'CTG',
	'M':'ATG', 'N':'AAC', 'O':'TAA', 'O':'TGA', 
	'P':'CCT', 'P':'CCA', 'P':'CCG', 'Q':'CAG',
  'R':'CGT', 'R':'CGC', 'R':'CGA', 'R':'CGG', 
  'R':'AGA', 'R':'AGG', 'S':'TCT', 'S':'TCC',
  'S':'TCA', 'S':'TCG', 'S':'AGT', 'S':'AGC',
  'T':'ACT', 'T':'ACC', 'T':'ACA', 'T':'ACG',
  'U':'CCC', 'V':'GTT', 'V':'GTC', 'V':'GTA',
  'W':'TGG', 'X':'GGG', 'Y':'TAT', 'Y':'TAC', 
  'Z':'CAA', '.':'TAG', ' ':'TTT', '#':'GTG'
}
'''

