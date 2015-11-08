Heuristic Search Algorithms
===========================

Overview
--------

Basically, there are three choices:

+ A*
+ IDA* (iterative deepening A*)
+ RBFS (recursive best first search)

If the state space is sufficiently small, use A*.
Otherwise, use IDA* or RBFS. 
If good solutions are spreaded among search pathes, use IDA*.
Otherwise, i.e., good solutions are condensed, use RBFS.

TODO
