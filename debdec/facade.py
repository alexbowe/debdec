import decode

class DeBruijnSequence:
  def __init__(self, window_size):
    # TODO: support larger alphabets as a parameter
    self.alphabet_size = 2 #len(alphabet) if alphabet else k
    self.window_size = window_size
    self.alphabet = range(self.alphabet_size)
    self._generator = None

  def __repr__(self):
    return "DeBruijnSequence(window_size=%s)" % (self.window_size)

  def __call__(self):
    if self._generator is None:
      self._generator = decode.make_generator(self.alphabet_size, self.window_size)
    return next(self._generator)

  def __iter__(self):
    return self

  def __len__(self):
    return self.alphabet_size**self.window_size

  def __contains__(self, s):
    if type(s) in {list, tuple}: return False
    return len(s) == self.window_size and len(set(s) & set(self.alphabet)) == self.alphabet_size

  def _int_idx(self, i):
    if not -len(self) <= i < len(self):
      raise IndexError("DeBruijnSequence index out of range")
    return decode.f(self.window_size, i)

  def _sequence_idx(self, s):
    #if s not in self: raise KeyError(s)
    return decode.decode_db(self.window_size, s)

  def _slice_idx(self, sl):
    start, stop, step = sl.indices(len(self))
    return [self._int_idx(i) for i in xrange(start, stop, step)]

  def __getitem__(self, key):
    dispatch = {int: self._int_idx,
                list: self._sequence_idx,
                tuple: self._sequence_idx,
                slice: self._slice_idx}[type(key)]
    return dispatch(key)

