from debdec import DeBruijnSequence

def wrapped_access(a, i, n):
  a = list(iter(a))
  return (a[:]+a[:n])[i:i+n]

def decode_test():
  for n in xrange(1, 9):
    yield (check_all_window_indexes, DeBruijnSequence(n))

# Goes through ALL indexes and checks
def check_all_window_indexes(dbs):
  for i in xrange(0, len(dbs)):
    assert check_window_index_invertibility(dbs, i)

# Gets block from dbs, then uses decode... this is not a unit test, but more like integration of two functions
# TODO: write a test that tests access() for a known sequence (k=5)
# TODO: write a test that tests generation for a known sequence
# TODO: write a test that tests decoding for a known sequence (with window extraction code tested as well)
# TODO: write randomised tests for higher dimensionality
def check_window_index_invertibility(dbs, idx):
  n = dbs.window_size
  block = wrapped_access(dbs, idx, n)
  result = dbs.decode(block)
  print result, idx, block
  if result != idx:
    print "Block: %s, expected index: %d, decoded index: %d" % (block, idx, result)
  return result == idx
