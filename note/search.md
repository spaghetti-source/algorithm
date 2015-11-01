Search Algorithms
=================

Overview
--------

1) A*
2) IDA* (iterative deepening A*)
3) RBFS (recursive best first search)

If the state space is sufficiently small, use A*.
Otherwise, use IDA* or RBFS.
If good solutions are spreaded among search pathes, use IDA*.
Otherwise, i.e., good solutions are condensed, use RBFS.

TODO
