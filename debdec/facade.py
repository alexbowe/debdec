import decode

# TODO: add documentation

class DeBruijnSequence:
  def __init__(self, window_size, alphabet_size=2):
    # TODO: support larger alphabets as a parameter
    self.alphabet_size = alphabet_size
    self.window_size = window_size
    self.alphabet = range(self.alphabet_size)

  def __repr__(self):
    c_name = self.__class__.__name__
    return "%s(window_size=%d, alphabet_size=%d)" % (c_name, self.window_size, self.alphabet_size)

  def __iter__(self):
    gen = decode.make_generator(self.alphabet_size, self.window_size)
    from itertools import islice
    return islice(gen, len(self))

  def __str__(self):
    return str(list(iter(self)))

  # They are defined as cyclic, but for ease of use I defined
  # fixed length. just add the period (length) if you need all occurences or something...
  # or use itertools cycle if you need to cycle it...
  def __len__(self):
    return self.alphabet_size**self.window_size

  def _int_idx(self, i):
    # TODO: implement this for alphabet size other than 2 and 4
    if self.alphabet_size != 2: raise NotImplementedError("Not yet implemented for alphabet size other than 2.")
    if not -len(self) <= i < len(self):
      raise IndexError("DeBruijnSequence index out of range")
    # TODO: test this... seems to not quite work
    return decode.f(self.window_size, i)

  # doesnt work like find(). If people want to find non-kmer substrings, they can implement it themselves
  def decode(self, s):
    if len(set(s) - set(self.alphabet)) != 0 or len(s) != self.window_size:
      return KeyError
    if self.alphabet_size == 2:
      return decode.decode_db(self.window_size, s)
    elif self.alphabet_size == 4:
      return decode.make_decoder(self.window_size)(s)
    # TODO: implement this for alphabet size other than 2 and 4
    raise NotImplementedError("Not yet implemented for alphabet size other than 2 and 4.") 

  def _slice_idx(self, sl):
    start, stop, step = sl.indices(len(self))
    return [self._int_idx(i) for i in xrange(start, stop, step)]

  # The semantics could change, so I'm considering DeBruijnSequence to NOT be a set of kmers
  # but just a sequence (as the name implies), so I don't use this function for reverse lookup (like a dictionary)
  # of a kmer.
  def __getitem__(self, key):
    dispatch = {int: self._int_idx,
                slice: self._slice_idx}[type(key)]
    return dispatch(key)

