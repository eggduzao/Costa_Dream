#!/usr/bin/env zsh

# Export
export PATH=$PATH:/hpcwork/izkf/bin/:/hpcwork/izkf/opt/meme_4.11.1/bin/:/hpcwork/izkf/opt/perl-5.24.0/bin/:
export PYTHONPATH=$PYTHONPATH:/hpcwork/izkf/lib64/python2.7/:/hpcwork/izkf/lib64/python2.7/site-packages/:/hpcwork/izkf/lib/python2.7/:/hpcwork/izkf/lib/python2.7/site-packages/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/hpcwork/izkf/lib64/:/hpcwork/izkf/lib/:/usr/lib64/R/lib:
export PERL5LIB=/hpcwork/izkf/opt/perl-5.24.0/lib/5.24.0/

python /work/eg474423/Costa_Dream/exp/meme_chip_positive/createMemeChIPPositive.py $1 $2 $3 $4 $5 $6 $7 $8


