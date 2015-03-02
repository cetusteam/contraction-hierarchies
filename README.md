Contraction Hierarchies
=======================
Version 2008/06/24 Revision 447

Please contact me (robert-geisberger@web.de) if you are using the source code
as I am very interested in the possible applications that arise from using it.

Introduction
------------

This program allows to create and test contraction hierarchies (CH), a
hierarchical, shortcut-based speedup technique for fast and simple
routing in road networks. The code was created as part of a diploma
thesis at Universitaet Karlsruhe (TH), Institute for Algorithmics II.

Papers: - R. Geisberger, P. Sanders, D. Schultes, D. Delling.
          Contraction Hierarchies: Faster and Simpler Hierarchical Routing
          in Road Networks. In 7th International Workshop on Efficient and
          Experimental Algorithms (WEA), 2008.

Required libraries
------------------

- STL
- Boost IO-streams & regular expressions (only for input/output)

How to run
----------

    cmake .
    make
    ./bin/node-order -s -p -f docu/example.ddsg -o output1-1.hcn -E output1-2.dsg -Z output1-3.sgr
    ./many/manyToMany output1-3.sgr 5 5 255


Sources of documentation
------------------------

In docs/index.html, there are some examples what command-line arguments
can be used to create a CH and perform tests.

Contact
-------
Robert Geisberger (geisberger@ira.uka.de)

Additional programs
-------------------

In the subfolder many/ is a test program for fast many-to-many distance table
computations. However, the source code and the command-line arguments are well
not documented well. It is an implementation of the algorithm desribed in the
paper

S. Knopp, P. Sanders, D. Schultes, F. Schulz, and D. Wagner. Computing
many-to-many shortest paths using highway hierarchies. In Workshop on
Algorithm Engineering and Experiments (ALENEX), 2007.

but using contraction hierarchies.
