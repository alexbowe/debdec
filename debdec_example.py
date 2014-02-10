from string import maketrans

import pyjudy

from debdec import make_decoder, sliding_windows

def reverse_complement(s):
  return [[3,2,1,0][x] for x in s[::-1]]

if __name__ == '__main__':
  import sys

  if len(sys.argv) != 2:
    print "usage: %s <k-mer length>" % (sys.argv[0])
    exit()
  try:
    k = int(sys.argv[1])
  except:
    print "usage: %s <k-mer length>" % (sys.argv[0])
    exit()
  if k < 1:
    print "k-mer length must be at least 1."

  p = True

  d = make_decoder(k)

  table = {base:idx for idx,base in enumerate("ACGT")}
  table_inv = [c for c,i in sorted(table.items(), key=lambda x:x[1])]

  j = pyjudy.Judy1Int()

  for line in sys.stdin:
    i = [table[c] for c in line.rstrip()]
    kmers = sliding_windows(i, k)
    for kmer in kmers:
      j.Set(d(kmer))
      j.Set(d(reverse_complement(kmer)))

  print list(j.iterkeys())
