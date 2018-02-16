Compile and Run on HCC (Crane)
------------------------------
Compile:

		module load compiler/gcc/7.1

		g++ -o kmerEst kmerCountEstimate.cpp -std=c++11 -O3 -march=native

Run:

		./kmerEst in.fasta kmer-Len sample-size 
  
Output:
  
  
  
