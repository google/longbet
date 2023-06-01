#! /bin/bash
echo Building R
cd ../../
R CMD REMOVE longBet
R CMD INSTALL longBet
cd longbet/tests/
echo Testing R
Rscript test_longBet.R
