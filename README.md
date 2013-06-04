db-decode
=========

An implementation of decodable [de Bruijn sequences][debseq]. Implementation of binary
de Bruijn decoding based on Knuth 4A[^knuth4A], whereas higher alphabets follow the
extended work by ...

As far as I can tell, this is the first public implementation of these algorithms.


De Bruijn Sequences
-------------------

Add definition and example.


Decodable De Bruijn Sequences
-----------------------------

Explain that if it is generated in a certain way, can be decoded quickly.
That is, given a n-length substring, we can calculate the position without
a complete lookup table.

Can also go from position to symbol.


Usage
-----

db-decode -n 26 -a 4 file
cat file | db-decode -n 26 -a 2


License
-------

Free to do whatever. All I ask is that if you make an improvement, contact me and submit
a pull request, and if you use this for a publication 
please put a link and my name somewhere.

[debseq]
