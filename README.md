debdec
======

**Version:** 0.1

An implementation of decodable [de Bruijn sequences][debseq]. Implementation of binary
de Bruijn decoding based on Knuth 4A, whereas higher alphabets follow the
extended work by Jonathan Tuliani.

As far as I can tell, this is the first public implementation of these algorithms, 
possibly even the first implementation.


Usage
-----

The included example script (which assumes an alphabet of size 4: {A,C,G,T} - currently only
binary and quaternary are supported) accepts one parameter for the subsequence lengths,
and accepts input from `stdin`. 

As an example, if `file` contained 26-length strings consisting of the symbols A, C, G and T:

    cat file | python debdec_example.py 26

If you want to use it as a library in your own code:

    >>> from debdec import make_decoder
    >>> d = make_decoder(5) # Assumes alphabet size of 4
    >>> d([0,1,0,1,1])
    352

I need to fix the interface to be cleaner, but if binary decoding is needed, `decode_db` can be used:

    >>> from debdec import decode_db
    >>> decode_db(5, [0,1,0,1,1]) # first parameter is the length
    9


Todo
----

1. Write unit tests
2. Refactor
3. Provide a facade for a nicer
4. Pack the sequence representation (currently uses lists of ints, could use bitvectors)
5. Remove the unnecessary use of Euclids algorithm
6. Profile and add memoization and tables where it could help.
7. Provide vectorized version, perhaps using CUDA
8. ???
9. Profit


License
-------

Free to do whatever you like. All I ask is that if you make an improvement, contact me and submit
a pull request. If you use it for something, attribution is certainly welcome and appreciated :)


Thanks
------

Jonathan Tuliani, Don Knuth, Chris Mitchell

[debseq]: https://en.wikipedia.org/wiki/De_Bruijn_sequence
