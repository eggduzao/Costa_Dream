#!/bin/zsh

# Eduardo's Libraries
export CC=gcc
export JAVA_TOOL_OPTIONS=-Xmx12000m
export FPICFLAGS=-fPIC
module load gcc &> /dev/null
module load python &> /dev/null
export PATH=$PATH:/home/eg474423/Installation/meme/bin # MEME SUITE
export PATH=$PATH:/home/eg474423/Installation/HINT/bin # HINT
#export PATH=$PATH:/home/eg474423/Installation/motifanalysis/bin # Motif Analysis
export PATH=$PATH:/home/eg474423/Installation/wiggle_to_psp_utilities # Cuellar Method
export PATH=$PATH:/home/eg474423/Installation/cutadapt-1.10/bin # cutadapt
export PATH=$PATH:/home/eg474423/Installation/FastQC # fastqc
export PATH=$PATH:/home/eg474423/Installation/trim_galore_zip/ # trim galore
export PATH=$PATH:/hpcwork/izkf/bin
export PATH=$PATH:/hpcwork/izkf/opt/bin
export PATH=$PATH:/home/eg474423/perl/bin
export PATH=$PATH:/work/eg474423/eg474423_Projects/trunk/TfbsPrediction/Code/bin
export PATH=$PATH:/hpcwork/izkf/backup/bin
export PYTHONPATH=$PYTHONPATH:/home/eg474423/Installation/scikit-learn-0.17.1/lib/python2.7/site-packages # scikit-learn-0.17.1
#export PYTHONPATH=$PYTHONPATH:/home/eg474423/Installation/scikit-learn-0.15.2/lib/python2.7/site-packages # scikit-learn-0.15.2
#export PYTHONPATH=$PYTHONPATH:/home/eg474423/Installation/pysam-0.9.0/lib/python2.7/site-packages # pysam-0.9.0
export PYTHONPATH=$PYTHONPATH:/home/eg474423/Installation/pysam-0.8.3/lib/python2.7/site-packages # Pysam 0.8.3
export PYTHONPATH=$PYTHONPATH:/home/eg474423/Installation/hmmlearn-0.2.0/lib/python2.7/site-packages # hmmlearn-0.2.0
export PYTHONPATH=$PYTHONPATH:/home/eg474423/Installation/HINT/lib/python2.7/site-packages # HINT
#export PYTHONPATH=$PYTHONPATH:/home/eg474423/Installation/motifanalysis/lib/python2.7/site-packages # Motif Analysis
export PYTHONPATH=$PYTHONPATH:/home/eg474423/Installation/bx-python-0.7.3/lib/python2.7/site-packages # bx-python-0.7.3
export PYTHONPATH=$PYTHONPATH:/home/eg474423/Installation/pyDNase-0.1.7/lib/python2.7/site-packages/ # PyDNase 0.1.7
export PYTHONPATH=$PYTHONPATH:/home/eg474423/Installation/cutadapt-1.10/lib/python2.7/site-packages/ # cutadapt
export PYTHONPATH=$PYTHONPATH:/hpcwork/izkf/lib64/python2.7/site-packages
export PYTHONPATH=$PYTHONPATH:/hpcwork/izkf/lib/python2.7/site-packages
export PYTHONPATH=$PYTHONPATH:/hpcwork/izkf/lib_clusterold/lib64/python2.6/site-packages/
export PYTHONPATH=$PYTHONPATH:/hpcwork/izkf/lib_clusterold/lib/python2.6/site-packages/
export PYTHONPATH=$PYTHONPATH:/hpcwork/izkf/backup/lib/python2.7/site-packages/
export PYTHONPATH=$PYTHONPATH:/hpcwork/izkf/backup/lib64/python2.7/site-packages/
export CPATH=/hpcwork/izkf/include:$CPATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/hpcwork/izkf/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/hpcwork/izkf/lib64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/hpcwork/izkf/opt/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/hpcwork/izkf/opt/lib64
export R_LIBS_USER=$R_LIBS_USER:/home/eg474423/R
export PERL5LIB=$PERL5LIB:/home/eg474423/Installation/meme/lib/perl
export PERL5LIB=$PERL5LIB:/home/eg474423/perl/lib64/perl5
export PERL5LIB=$PERL5LIB:/home/eg474423/perl/lib/perl5
export PERL5LIB=$PERL5LIB:/home/eg474423/perl/share/perl5
export PERL5LIB=$PERL5LIB:/hpcwork/izkf/lib/perl

# Input
inLoc=$1
outLoc=$2
cell=$3

# Evaluate Norm Factor
bamFileName=$inLoc$cell"/"$cell"_DNase.bam"
hsFileName=$inLoc$cell"/DNase_Peaks.bed"
normName=$outLoc"norm/"$cell".txt"
python EvaluateNormFactor.py $bamFileName $hsFileName $normName

# Evaluate Dictionaries
k_nb="6"
allTagsFg="N"
bamFileName=$inLoc$cell"/"$cell"_DNase.bam"
hsFileName=$inLoc$cell"/DNase_Peaks.bed"
fastaFileName="/hpcwork/izkf/projects/TfbsPrediction/Data/HG19/hg19.fa"
csFileName="/hpcwork/izkf/projects/TfbsPrediction/Data/HG19/hg19.chrom.sizes.filtered"
dictName=$outLoc"dict/"$cell
python CreateKmerDict.py $k_nb $allTagsFg $bamFileName $hsFileName $fastaFileName $csFileName $dictName

# Evaluate Tables
python CreateKmerTable.py $k_nb "All" $dictName"_obsDict.p" $dictName"_expDict.p" $normName $outLoc"table/"$cell"_All.txt"
python CreateKmerTable.py $k_nb "F" $dictName"_obsDictF.p" $dictName"_expDictF.p" $normName $outLoc"table/"$cell"_F.txt"
python CreateKmerTable.py $k_nb "R" $dictName"_obsDictR.p" $dictName"_expDictR.p" $normName $outLoc"table/"$cell"_R.txt"

# Format Bias Tables
python FormatBiasTable.py


