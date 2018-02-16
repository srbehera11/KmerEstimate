compile on HCC:

module load compiler/gcc/7.1

g++ -o kmerEst kmerCountEstimate.cpp -std=c++11 

Run:

./kmerEst <in.fasta> <kmerLen>        
