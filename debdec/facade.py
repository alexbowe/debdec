import decode

class DeBruijnSequence:
  def __init__(self, n):
    self.k = 2#len(alphabet) if alphabet else k
    self.n = n
    self.alphabet = range(self.k)
    self._generator, self._decoder = None, None

  def __repr__(self):
    return "DeBruijnSequence(n=%s, k=%s)" % (self.n, self.k)

  def __call__(self):
    if self.generator is None:
      self.generator = decode
    return self.generator()

  def _int_idx(self, i):
    return decode.f(self.n, i)

  def _slice_idx(self, sl):
    return [self._int_idx(i) for i in xrange(sl.start, sl.stop, sl.step)]

  def _iter_idx(self, s):
    return decode.decode_db(self.n, s)

  def __getitem__(self, key):
    dispatch = {int: self._int_idx,
                list: self._iter_idx,
                tuple: self._iter_idx,
                slice: self._slice_idx}[type(key)]
    return dispatch(key)

  def __next__():
    pass

# array access for integers
# array access for strings
