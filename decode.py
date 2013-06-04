#! /usr/local/bin/python

import itertools as it
import math


## TABLES ##

db2      = [0,0,1,1]
sum_db2  = [0,0,0,1]
db2_tups = {(0,0):0, (0,1): 1, (1,1): 2, (1,0): 3}
db3      = [0,0,0,1,1,1,0,1]
sum_db3  = [0,0,0,0,1,0,1,1]
db3_tups = {(0,0,0):0, (0,0,1):1, (0,1,1):2, (1,1,1):3,
            (1,1,0):4, (1,0,1):5, (0,1,0):6, (1,0,0):7}
db4      = [0,0,0,0,1,0,1,0,0,1,1,1,1,0,1,1]
sum_db4  = [0,0,0,0,0,1,1,0,0,0,1,0,1,0,0,1]
db4_tups = {(0,0,0,0):0, (0,0,0,1):1, (0,0,1,0):2, (0,1,0,1):3,
            (1,0,1,0):4, (0,1,0,0):5, (1,0,0,1):6, (0,0,1,1):7,
            (0,1,1,1):8, (1,1,1,1):9, (1,1,1,0):10, (1,1,0,1):11,
            (1,0,1,1):12, (0,1,1,0):13, (1,1,0,0):14, (1,0,0,0):15}


## UTILITY FUNCTIONS ##

lg = lambda x: math.log(x, 2)

def windows(seq, length=2, overlap=0):
  seq=iter(seq)
  result = list(it.islice(seq, length))
  while len(result) == length:
    yield result
    result = result[length - overlap:]
    result.extend(it.islice(seq, length-overlap))

def sliding_windows(seq, length=2):
  return windows(seq, length, length - 1)

def rotate(x, i=1):
  x = list(x)
  return x[i:] + x[:i]

def rotate_left(x, i=1):
  return rotate(x, i)

def rotate_right(x, i=1):
  return rotate(x, -i)

def paired_xor(b):
  # if this was implemented using bitmaps, could be done with
  # rotate left, sub 1
  return [a ^ b for a,b in sliding_windows(b, 2)]

# Only works for the rightmost 32 bits
# can possibly be looped over blocks if support for longer values is required
# but in this context x is the length of a kmer (probably < 2^32)
def count_right_zero_bits(x):
  mod_bit_position = [ 32,  0,  1, 26,  2, 23, 27,  0,  3, 16,
                       24, 30, 28, 11,  0, 13,  4,  7, 17,  0,
                       25, 22, 31, 15, 29, 10, 12,  6,  0, 21,
                       14,  9,  5, 20,  8, 19, 18 ]
  return mod_bit_position[(-x & x) % 37]


## DECODING M-ARY DE BRUIJN SEQUENCES ##

# NOTE: could replace cycle with bit-twiddling
# and keep v as an integer, then not have to call psi_inv
def A_inv(v, t):
  v_prime = rotate_right(v, t)
  return psi_inv(v_prime)

def B_inv(v, t):
  n = len(v)
  b = A_inv(v,t)
  if sum(v[:-t]) == n - t:
    b = (b - 2**t) % 2**n
  return b

def C_inv(v, t):
  if t == 0:
    return psi_inv(v)
  n = len(v)
  left_w = v[n-t:]
  right_w = v[:-t]
  mid_w = [-sum(left_w + right_w) % 2]
  w = left_w + mid_w + right_w
  return psi_inv(w[-n:])

def D_inv(v, t):
  if t > 0:
    return D_inv(v, t-1)
  n = len(v)
  w_0 = [-sum(v) % 2]
  w = v + w_0
  return psi_inv(w[-n:])

def E_inv(v, t):
  return D_inv(v,t) if t == 0 else C_inv(v, t-1)

def F_inv(v, t):
  if t == 0 or t == 1:
    return E_inv(v, t)
  n = len(v)
  u = v[:-t+1]
  right_w = psi(n-t+1,(psi_inv(u) + 1)% 2**(n-t+1))
  left_w = v[-t+1:]
  mid_w = [-sum(left_w + right_w) % 2]
  w = left_w + mid_w + right_w
  return (psi_inv(w[-n:]) - 1) % (2**n)

# Maps a number in the ring Z_{2^n} to
# an n-bit sequence from Z_2 (for some n)
# (in other words, returns binary representation)
def psi(n, x):
  b = map(int, list(bin(x)[2:]))
  return [0] * (n - len(b)) + b

# Inverse of the above function: takes a n-bit sequence
# and convertes it back to an element of the composite alphabet
def psi_inv(b):
  return int(''.join(map(str,b)), 2)

def theta(n, x):
  code = psi(n, x)
  # add parity bit
  return [sum(code) % 2] + code

# to construct [s]
def d_x(n, k, l, x):
  return psi(n, x) * k + theta(n, x) * l

def decode_sdb(n, k, l, v, f):
  n = len(v)
  if l == 0:
    t = f % n
    return A_inv(v,t) if f < n*(k-1) else B_inv(v,t)
  elif k == 0:
    t = f % (n+1)
    return E_inv(v,t) if f < (n+1)*(l-1) else F_inv(v,t)
  elif f < n*k:
    t = f % n
    return A_inv(v,t) if f < n*(k-1) else C_inv(v,t)
  else:
    t = (f-n*k) % (n+1)
    return E_inv(v,t) if f <= n*k+(n+1)+(l-1) else D_inv(v,t)


## BINARY BASE DE BRUIJN DECODER (D + R) ##

# Position of 1^n
def j(n):
  if n <= 0: return
  if n <= 3: return n
  b = lambda t: (2**t ^ (-1)**t)/3
  t = count_right_zero_bits(n)
  q = n >> t
  if q == 1: return 2**(n-1) + b(t)
  bt = b(t + 2)
  m = int(lg(q))
  m += (4*q)/2**m == 5
  if m % 2 != 0: bt = -bt
  return 2**(n-1) + bt

def delta(n):
  assert n%2 == 0 # only for even n
  return 1 - 2 * (j(n/2) % 4 == 3)

# position of 1^(2n - 1)
def k_n(n):
  assert n >= 2
  if n == 2: return 6
  return (2**(n-1)-2)*(2**n+2) + 2*j(n) + 2

# position of (01)^(n-1)0
def l(n):
  assert n >= 2
  if n == 2: return 2
  jn = j(n)
  if jn % 4 == 1: return (2**(n-1)+1)*(jn - 1)
  elif jn % 4 == 3: return (2**(n-1)+1)*(2**n + jn - 3)
  assert False # undefined for jn % 4 == 2

def c(n, k):
  return f_pos(n, k//2) if k%2 == 0 else f_neg(n, (k-1)//2)

def f_pos(n, k):
  k = k % (2**n + 2)
  return f(n, k-1) if 0 < k <= j(n) + 1 else f(n, k-2)

def f_neg(n, k):
  k = k % (2**n - 2)
  return f(n, k+1) if 0 <= k < j(n) else f(n, k+2)

def sum_f_2n(n, k):
  if k == 0: return 0
  jn,ln,kn = j(n), l(n), k_n(n)
  if jn % 4 < 3:
    if        0 < k <= ln + 2:   return sum_c(n, k-1)
    elif ln + 2 < k <= kn + 3:   return 1 + sum_c(n, k-3)
    elif kn + 3 < k <= 2**(2*n): return sum_c(n, k-4)
  elif jn % 4 == 3:
    if        0 < k <= kn + 1:   return sum_c(n, k-1)
    elif kn + 1 < k <= ln + 3:   return 1 + sum_c(n, k-2)
    elif ln + 3 < k <= 2**(2*n): return sum_c(n, k - 4)

def sum_f_n_plus_1(n, k):
  if k == 0: return 0
  jn = j(n)
  if 0 < k <= jn or 2**n + jn < k <= 2**(n+1): return sum_sum_f_2n(n/2, k)
  elif jn < k <= 2**n + jn: return 1 + k + sum_sum_f_2n(n/2, k + delta(n))

def sum_f(n, k):
  k = k % 2**n
  if k == 0: return 0
  if n == 2: return sum_db2[k]
  if n == 3: return sum_db3[k]

  # f_2n
  if n%2 == 0:
    r = sum_f_2n(n/2, k)
    assert r != None
    return r
  else: # f_2n+1
    r = sum_f_n_plus_1(n - 1, k)
    assert r != None, "n = %d, k = %d" %(n, k)
    return r

def sum_c(n, k):
  return sum_f_pos(n, int(math.ceil(k/2.0))) + sum_f_neg(n, k//2)

def sum_sum_c(n, k):
  return sum_f_pos(n, k//2) if k % 2 == 0 else sum_f_neg(n, (k-1)//2)

def sum_sum_f_2n(n, k):
  if k == 0: return 0
  jn, ln, kn = j(n), l(n), k_n(n)
  offset = (jn % 4 < 3) * k//2**(2*n)
  k = k % 2**(2*n)
  if k == 0: return offset
  if jn % 4 < 3:
    if      0 < k <= ln + 2:   return offset + sum_sum_c(n, k-1)
    if ln + 2 < k <= kn + 3:   return offset + 1 + k + sum_sum_c(n, k-3)
    if kn + 3 < k <= 2**(2*n): return offset + sum_sum_c(n, k-4)
  if jn % 4 == 3:
    if      0 < k <= kn + 1:   return offset + sum_sum_c(n, k-1)
    if kn + 1 < k <= ln + 3:   return offset + 1 + k + sum_sum_c(n, k-2)
    if ln + 3 < k <= 2**(2*n): return offset + 1 + sum_sum_c(n, k-4)

def sum_f_pos(n, k):
  offset = k//(2**n + 2)
  k = k % (2**n + 2)
  return offset + (k > 0) * (sum_f(n, k-1) if 0 < k <= j(n) + 1 else 1 + sum_f(n, k-2))

def sum_f_neg(n, k):
  offset = k//(2**n - 2)
  k = k % (2**n -  2)
  return offset + (sum_f(n, k+1) if 0 <= k < j(n) else 1 + sum_f(n, k+2))

def f_n_plus_1(n, k):
  assert n%2 == 0
  jn = j(n)
  if 0 < k <= jn or 2**n + jn < k <= 2**(n+1):  return sum_f(n, k)
  elif jn < k <= 2**n + jn: return 1 + sum_f(n, k + delta(n))

def f_2n(n, k):
  jn, ln, kn = j(n), l(n), k_n(n)
  if jn % 4 < 3:
    if        0 < k <= ln + 2:   return c(n, k-1)
    elif ln + 2 < k <= kn + 3:   return c(n, k-3)
    elif kn + 3 < k <= 2**(2*n): return c(n, k-4)
  elif jn % 4 == 3:
    if        0 < k <= kn + 1:   return c(n, k-1)
    elif kn + 1 < k <= ln + 3:   return c(n, k-2)
    elif ln + 3 < k <= 2**(2*n): return c(n, k-4)

# Access element of db
def f(n, k):
  assert n >= 1
  if n == 1: return k % 2
  k = k % 2**n
  if k < n: return 0
  if n == 2: return db2[k]
  if n == 3: return db3[k]

  if n % 2 == 0:
    return f_2n(n/2, k) % 2
  else:
    return f_n_plus_1(n-1, k) % 2

def decode_db(n, u):
  assert n >= 1
  if n == 1: return u[0]
  if n == 2: return db2_tups[tuple(u)]
  if n == 3: return db3_tups[tuple(u)]
  if n == 4: return db4_tups[tuple(u)]

  if n % 2 == 1: # -> n = m + 1
    m = n - 1
    p = decode_db(m, paired_xor(u))
    d_m, j_m = delta(m), j(m)
    k = p - d_m * (p >= j_m) + 2**m * (p == j_m) * (d_m == 1)
    if f(n, k) == u[0]: return k
    elif p < j_m:  return k + 2**m - d_m
    elif p == j_m: return k + d_m
    elif p > j_m:  return k + 2**m + d_m
  else: # -> n = 2m
    m = n//2
    pos = lambda k: k + 1 + (k > j(m))
    neg = lambda k: k - 1 - (k > j(m))
    j2n, ln, kn = j(n), l(m), k_n(m)
    if u == [0] * n: return 0
    if u == [1] * n: return j2n
    if u == [0,1] * (m): return ln + 1 + (kn < ln)
    if u == [1,0] * (m): return ln + 2 + (kn < ln)
    alpha = [a for i,a in enumerate(u) if i % 2 == 0]
    beta = [a for i,a in enumerate(u) if i % 2 == 1]
    k0, k1 = decode_db(m, alpha), decode_db(m, beta)
    k0_pos, k1_pos = pos(k0), pos(k1)
    if alpha == [0]*m:
      (x,y) = isolve(2**m+2,2-2**m,neg(k1) - k0_pos + ((pos(k0) - pos(k1)) % 2 != 0))
      p = (2*(neg(k1) + y*(2**m - 2))) % (2**n-4)
    elif beta == [0] * m:
      (x,y) = isolve(2**m-2,-2-2**m,k1_pos - neg(k0) -1 - ((pos(k0) - pos(k1)) % 2 == 0))
      p = (2*(neg(k0) + x*(2**m - 2))+1) % (2**n-4)
    elif alpha == [1] * m:
      (x,y) = isolve(2**m+2,2-2**m,neg(k1) - k0_pos - ((pos(k0) - pos(k1)) % 2 != 0))
      p = (2*(neg(k1) + y*(2**m - 2))) % (2**n-4)
    elif beta == [1]*m:
      (x,y) = isolve(2**m-2,-2-2**m,k1_pos - neg(k0) -1 + ((pos(k0) - pos(k1)) % 2 == 0))
      p = (2*(neg(k0) + x*(2**m - 2))+1) % (2**n-4)
    elif (pos(k0) - pos(k1)) % 2 == 0:
      (x,y) = isolve(2**m+2,2-2**m,neg(k1)-k0_pos)
      p = (2*(k0_pos + x*(2**m + 2))) % (2**n-4)
    else:
      (x,y) = isolve(2**m-2,-2**m-2,k1_pos - neg(k0) - 1)
      p = (2*(neg(k0) + x*(2**m - 2)) + 1) % (2**n-4)
    return (p + 1 + (p >= kn) + 2*(p >= ln) )

# www.math.utah.edu/~carlson/hsp2004/PythonShortCourse.pdf
# uses a variation of Euclid's algorithm
def isolve(a, b, c):
  q, r = divmod(a, b)
  if r == 0:
    return (0, c/b)
  else:
    u, v = isolve(b, r, c)
    return (v, u - q*v)

def find_k_and_l(n):
  # solve diophantine eqn for 2^n = nk + (n+1)l
  # n >= 1, k,l >= 0
  a, b, c = n, n+1, 2**n
  x, y = isolve(a, b, c)
  # find positive/zero values
  t = 0 if x >= 0 and y >= 0 else 1 + (-x - 1)/b
  k, l = x + t*b, y - t*a
  assert k != 0 or l != 0
  assert k >= 0 and l >= 0
  assert a*k + b*l == c
  # might want to yield these later instead, and increase t until the
  # conditions above aren't met
  return (x + t*b, y - t*a)

# assumes c = 4, y = z = 2
def make_decoder(n, k=None, l=None):
  if k is None or l is None:
    k,l = find_k_and_l(n)

  assert 2**n == n*k + (n+1)*l

  def decode(w):
    assert len(w) == n
    u = [psi(2,x)[0] for x in w]
    v = [psi(2,x)[1] for x in w]
    f = decode_db(n, u)
    g = decode_sdb(n, k, l, v, f)
    return g*2**n + f

  return decode


## DE BRUIJN SEQUENCE GENERATORS ##

# Trivial de Bruijn cycle for a window size of one (n = 1)
def one_cycle(m = 2):
  return it.cycle(xrange(m))

# Knuth Vol 4A pg 302 (54): m-ary de Bruijn cycle with n=2
def two_cycle(m = 2):
  while True:
    for i in xrange(m):
      for x in xrange(i):
        yield x
        yield i
      yield i

# Knuth's suggestion for r (alg. D) such that
# 0 <= r <= m and gcd(m^n - r, m^n + r) = 2
# (r = 1 when m is odd, r = 2 when m is even)
def choose_r(m):
  return (m + 1) % 2 + 1

# Algorithm D (Doubly recursive de Bruijn cycle generation)
# Page 304 Knuth Volume 4A
def algorithm_D(f, f_prime, m, n, r=None):
  r = r or choose_r(m)
  # D0 [Initialize variables.]
  x = x_prime = m
  t = t_prime = y = y_prime = 0
  # D1 [Possibly invoke f.]
  while True:
    if t != n or x >= r: y = f.next()
    # D2 [Count repeats.]
    if x != y:
      x = y
      t = 1
    else:
      t = t + 1
    # D3 [Output from f.]
    while True:
      yield x
      # D4 [Invoke f'.]
      while True:
        y_prime = f_prime.next()
        # D5 [Count repeats.]
        if x_prime != y_prime:
          x_prime = y_prime
          t_prime = 1
        else:
          t_prime = t_prime + 1
        # D6 [Possibly reject f'.]
        if not (t_prime == n and x_prime < r and (t < n or x_prime < x)):
          break
      if t_prime == n and x_prime < r and x_prime == x:
        continue # goto D3
      # D7 [ Output from f'.]
      yield x_prime
      if not (t_prime == n and x_prime < r):
        break

def algorithm_R(f, m, n):
  x = y = t = 0
  # R1 [Output.]
  while True:
    yield x
    # R2 [Invoke f.]
    if x == 0 or t < n:
      y = f.next()
    while True:
      # R3 [Count ones]
      t = t + 1 if y == 1 else 0
      # R4 [Skip one?]
      if t == n and x != 0: y = f.next()
      else: break
    # R5 [Adjust x.]
    x = (x + y) % m

if __name__ == '__main__':
  import sys

  if len(sys.argv) != 2:
    print "usage: %s <k-mer length>"
  try:
    k = int(sys.argv[1])
  except:
    print "usage: %s <k-mer length>"
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
