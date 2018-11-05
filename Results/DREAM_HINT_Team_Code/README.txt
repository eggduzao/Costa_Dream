
######################################################################################
# ENCODE-DREAM in vivo transcription factor binding site (TFBS) Prediction Challenge
# Team: HINT
# Authors: Eduardo G. Gusmao, Zhijian Li and Ivan G. Costa
# Email: eduardo.gusmao@rwth-aachen.de / ivan.costa@rwth-aachen.de
######################################################################################

######################################################################################
# 1. Prerequisites
######################################################################################

To run our solution to the challenge, please make sure you have the following tools and python packages:

	1.1. numpy (http://www.numpy.org/) [tested in version 1.10.4]
	1.2. bedtools (http://bedtools.readthedocs.io/en/latest/) [tested in version 2.25.0]
	1.3. BioPython (http://biopython.org/) [tested in version 1.66]
	1.4. pysam (https://github.com/pysam-developers/pysam) [tested in version 0.8.3]
	1.5. MOODS (http://www.regulatory-genomics.org/wp-content/uploads/2016/04/MOODS_1.0.1.tar.gz)

The MOODS version needs to be the 1.0.1 as in the link provided above. To install, unzip the file and simply type the following commands:

cd <MOODS_DIRECTORY>/src
make
cd <MOODS_DIRECTORY>/python
sudo python setup.py install

######################################################################################
# 2. Execution Instructions:
######################################################################################

To execute our solution script, please open the file "main.py" in any text editor. Modify the lines marked with the comment "Change this line" at the end.

	2.1. The variable "originalDNaseLocation" should point to the path in which all training DNase-seq datasets reside.
	2.2. The variable "genomeFileName" should point to the fasta file containing the genome sequence.
	2.3. The variable "annotationBedFileName" should point to the bed file containing the non-merged challenge region annotations

After changing these lines just execute the following command:

python main.py

The output for the specified cell type (variable "cell" in the code) and transcription factor (variable "factor" in the code) will be gzipped in the correct format in the folder "out"


