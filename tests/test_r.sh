#! /bin/bash
echo Building R
cd ../../
R CMD REMOVE longbet
R CMD INSTALL longbet
cd longbet/tests/simulation
echo Testing R
Rscript test_dgp.R
