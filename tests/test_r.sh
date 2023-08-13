#! /bin/bash
echo Building R
cd ../../
R CMD REMOVE longbet
R CMD INSTALL longbet
cd longbet/tests/
echo Testing R
Rscript test_longbet.R
