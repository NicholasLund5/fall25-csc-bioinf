# Deliverable 3


Step 1 (Python setup):
   Get python/cython files files, remove unessesary data, convert cython to python, change data/import locations, add a timer for tests, add pytest to the github CI install

Step 2 (Codon Port):
   Following very similar process to deliverables 1 and 2
   The tree class in tree.codon gave me a pretty hard time with compilation errors. Codon's strict type checking was a pain point here.

Step 3 (Evalutate.sh):
   Very similar to deliverable 1, just run both python and codon tests


All tests:
test_distances
test_upgma
test_neighbor_joining
Have been ported and are working in both Codon and Python

Time to complete: 8 hours (ish)