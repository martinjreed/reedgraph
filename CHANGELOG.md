#Change log
## [1.48] -2014-09-21
### Changes
- rg.max.concurrent.flow.int now uses prescaling
### Fixed
- Fixed returning lambda from min_congestion flows
- Commented out unused variables

## [1.47] - 2014-04-10 
### Changes

- Now includes reedgraphEnv to set BLOOMLENGTH

## [1.44] - 2014-04-09
- Minor bug fixes

- max\_concurrent\_flow\_int optimises for min average gamma if gamma equal

## [1.42]

- Two nasty bugs in max-flow:max\_concurrent\\_flow\_int fixed

## [1.4]
- now rg.max.concurrent.flow.int finds best gamma

- Minor changes in max-flow.cpp (copying map in max\_concurrent\_flow\_int line 615)

- Version 1.3 Integer max flow tested (see rg.test.int.versus.nonint.flow), may not compile on Windows

- Updated max-flow.cpp so that it compiles

## [1.1]
- working but some bugs

- * R/reedgraph-algorithms.R: Added new, cleaner boost version of max commodity flow

- Updated all to create a package. All works on Mac and WIndows

- Added license (GPL) and split igraph code out WIll be removing igraph code in next version

- * R/reedgraph-algorithms.R: update to separate igraph out Should work cross-platform now
Will remove igraph after this as it is too slow (and buggy)

- * R/reedgraph-algorithms.R: small changes

- * R/reedgraph-bloom.R: New file to do Bloom Filter forwarding

- Updated updateExplicitFlow and elsewhere to cope with node names other than the same as the indices of the nodes.
Added st-cut (using approximate single flow)

- Code updated for PURSUIT 2012 ICC submission

- Updated for changes in Rcpp after R 2.13

- Added eInternal and extra testing between phases for lowest gamma using randomly selected path (between current and last phase)

- Tidied up code, an changed to "lowest" in analyse.inner

- Few bug fixes (fixed problem with permutation="fixed")
New function in reedgraph-scratch.R : rg.try.single.demands probably
does not work. Need to remove?

- Bug fixes for permutation added 

- Added permutation to demand order
Now can select random, lowest (by path cost), natural order or specify
a vector for the order.


- This version of SEXP rg\_max\_concurrent\_flow\_int\_c() produces very promising results. It found all of the best gamma solutions for a test set of more than 100 cases (albeit small cases that were also brute force tested to make sure it found the lowest).
Still needs more investigation though.

- This is more or less stable version, includes code in rg\_fleischer\_max\_concurrent\_flow\_restricted\_c for searching for best gamma. This cannot work with a single length vector! Hence, removing to fork the rg\_fleischer\_max\_concurrent\_flow\_restricted(\_c) code with a different approach.

- Removed demand.edge\_flow\_map as it is not needed.

- This version has rg.max.concurrent.flow complete however, it does not use a sensible value for scaled demands, so Beta may be less than 1 (not valid for the algorithm) or very large (very slow).

