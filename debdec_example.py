from debdec import make_decoder, sliding_windows

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

  for line in sys.stdin:
    i = [table[c] for c in line.rstrip()]
    kmers = sliding_windows(i, k)
    for kmer in kmers:
      print d(kmer), "".join([table_inv[i] for i in kmer]) if p else ""
