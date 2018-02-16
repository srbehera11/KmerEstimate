Compile on HCC:

		module load compiler/gcc/7.1

		g++ -o kmerCountEstimate kmerCountEstimate.cpp -std=c++11 -fopenmp -O3 -march=native

Run:

		./kmerEst in.fasta kmer-Len sample-size 
  
Output:
  
  
  
