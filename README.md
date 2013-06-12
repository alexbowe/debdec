debdec
======

An implementation of decodable [de Bruijn sequences][debseq]. Implementation of binary
de Bruijn decoding based on Knuth 4A[^knuth4A], whereas higher alphabets follow the
extended work by Jonathan Tuliani[^Tuliani].

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

as lib:

arity defaults to 2

make_decoder(arity=4, order=25)
make_coroutine(arity=4, order=25) to generate sequence in order O(1) time?
make_encoder(arity, order)


License
-------

Free to do whatever. All I ask is that if you make an improvement, contact me and submit
a pull request, and if you use this for a publication 
please put a link and my name somewhere.


Thanks
------

Tuliani, Knuth, Chris...

[debseq]: https://en.wikipedia.org/wiki/De_Bruijn_sequence
