//#include <google/sparse_hash_map>
//#include <google/dense_hash_map>
#include "MurmurHash3.cpp"
#include <iostream>
#include <climits>
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <time.h>
#include "metrohash64.cpp"
#include <stdint.h>
#include <unordered_map>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <cassert>
#include <string.h>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <math.h>
#include <sys/time.h>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <list>
#include <stack>
#include <limits.h>
#include <map>
#include <bitset>
#include <ctime>
#include <queue>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <cstring>
#include <iostream>
#include <random>
#include <cinttypes>
//#include "dna_test.h"
#include "ntHashIterator.hpp"

#define SPP_MIX_HASH 1
#include "sparsepp/spp.h"

using spp::sparse_hash_map;

using namespace std;
//KSEQ_INIT(gzFile, gzread)
KSEQ_INIT(int, read)
std::map<char, char> mapp = {{'A', 'T'}, {'C', 'G'}, {'G', 'C'}, {'T', 'A'}, {'N', 'N'}};

// Function to reverse string and return reverse string pointer of that
/*void ReverseConstString(char *str)
{
    int start, end, len;
    char temp, *ptr = NULL;
    len = strlen(str);  
    ptr = (char*) malloc(sizeof(char)*(len+1)); 
    strcpy(ptr,str);           
    for (start=0,end=len-1; start<=end; start++,end--)
    {
        temp = ptr[start];
        ptr[start] = mapp.at(ptr[end]);       
        ptr[end] = mapp.at(temp);
    }
    ptr[len] = '\0';
    strcpy(str, ptr);
    free(ptr); 
}*/

/*char complement(char n)
{   
    switch(n)
    {   
    case 'A':
        return 'T';
    case 'T':
        return 'A';
    case 'G':
        return 'C';
    case 'C':
        return 'G';
    }   
    assert(false);
    return ' ';
}*/   
static const int MultiplyDeBruijnBitPosition[32] =
{
  0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
  31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
};
unsigned trailing_zeros(unsigned n) {
    return n ? __builtin_ctz(n) : -1;
}

static const char basemap[255] =
    {   
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*   0 -   9 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  10 -  19 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  20 -  29 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  30 -  39 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  40 -  49 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  50 -  59 */
        '\0', '\0', '\0', '\0', '\0',  'T', '\0',  'G', '\0', '\0', /*  60 -  69 */
        '\0',  'C', '\0', '\0', '\0', '\0', '\0', '\0',  'N', '\0', /*  70 -  79 */
        '\0', '\0', '\0', '\0',  'A',  'A', '\0', '\0', '\0', '\0', /*  80 -  89 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0',  't', '\0',  'g', /*  90 -  99 */
        '\0', '\0', '\0',  'c', '\0', '\0', '\0', '\0', '\0', '\0', /* 100 - 109 */
        '\0', '\0', '\0', '\0', '\0', '\0',  'a',  'a', '\0', '\0', /* 110 - 119 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 120 - 129 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 130 - 139 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 140 - 149 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 150 - 159 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 160 - 169 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 170 - 179 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 180 - 189 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 190 - 199 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 200 - 209 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 210 - 219 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 220 - 229 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 230 - 239 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 240 - 249 */
        '\0', '\0', '\0', '\0', '\0'                                /* 250 - 254 */
    };

/*uint32_t murmur3_32(const uint8_t* key, size_t len, uint32_t seed) {
  uint32_t h = seed;
  if (len > 3) {
    const uint32_t* key_x4 = (const uint32_t*) key;
    size_t i = len >> 2;
    do {
      uint32_t k = *key_x4++;
      k *= 0xcc9e2d51;
      k = (k << 15) | (k >> 17);
      k *= 0x1b873593;
      h ^= k;
      h = (h << 13) | (h >> 19);
      h = (h * 5) + 0xe6546b64;
    } while (--i);
    key = (const uint8_t*) key_x4;
  }
  if (len & 3) {
    size_t i = len & 3;
    uint32_t k = 0;
    key = &key[i - 1];
    do {
      k <<= 8;
      k |= *key--;
    } while (--i);
    k *= 0xcc9e2d51;
    k = (k << 15) | (k >> 17);
    k *= 0x1b873593;
    h ^= k;
  }
  h ^= len;
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;
  return h;
}
*/
/*
struct comparator {
 bool operator()(uint64_t i, uint64_t j) {
 return i > j;
 }
};

struct CompareBySecond {
    bool operator()(pair<int, int> const & a,
                              pair<int, int> const & b)
    { return a.second > b.second; }
};
*/
unsigned trailing_zeros(uint64_t n) {
    return n ? __builtin_ctzll(n) : -1;
}

void printHelp()
{
    
    cout << "KmerEst [options] -f <fasta/fastq> -k <k-mer length>  -s <sample size> -o <output file>"    << endl
    << "  -h               help"                                   << endl
    << "  -f <file>       Input sequence file "                << endl
    << "  -k <k-mer size >        kmer size (default 31) "        << endl
    << "  -s <sample size>        sample size (default 25m)"        << endl
     << "  -c coverage>       coverage (default 64)"        << endl
    << "  -o         	  Prefix of the Output file " << endl;
    
    exit(0);
}
 
int main(int argc, char** argv)
{
    
    if(argc == 1){
      cout << argv[0] << " -f <seq.fa> -k  <kmerLen> -s <minHeap_Size> -c <coverage> -o <out.txt>" << endl;
      exit(0);
    } 
    int n=31;
    int s = 25000000;
    int cov = 64;
    string f = "", outf = "";
    for (int c = 1; c < argc; c++)
        {
            
            if (!strcmp(argv[c], "-h"))       { printHelp(); }
            else if (!strcmp(argv[c], "-k"))     { n = atoi(argv[c+1]); c++; }
            else if (!strcmp(argv[c], "-f"))    { f = argv[c+1]; c++; }
            else if (!strcmp(argv[c], "-s"))    { s = atoi(argv[c+1]); c++; }
            else if (!strcmp(argv[c], "-c"))    { cov = atoi(argv[c+1]); c++; }
            else if (!strcmp(argv[c], "-o")) { outf = argv[c+1]; c++; }
        }
        
       if (f.empty()  || outf.empty())
        {
          printHelp();
        }
    //int n = atoi(argv[2]);
    FILE *fp;
    kseq_t *seq;
    int l;
    //fp = fopen(argv[1], "r");
    fp = fopen(f.c_str(), "r");
    if( fp == Z_NULL){
      cout<<"File: "<< f  << " does not exist" <<endl;
      exit(1);
    }
    seq = kseq_init(fileno(fp));
    int k = s;
    // int k = atoi(argv[3]); // size of maxheap i.e. sample size
    //unordered_map<uint64_t, pair<uint8_t, uint32_t>> MAP; // (<hash>, <tz, count>>) 
    //sparse_hash_map<uint32_t, pair<uint8_t, uint32_t>> MAP;
    //vector <google::sparse_hash_map<uint32_t, pair<uint8_t, uint32_t>> > MAP(64);
    //MAP.set_empty_key(-1);
    // for(int i=0; i<64; i++) MAP[i].set_deleted_key(0);
    typedef sparse_hash_map<uint64_t, uint32_t> SMap;
    vector<SMap> MAP(64);
   // vector < sparse_hash_map<uint32_t, uint32_t> > MAP(64);
    cout << "read the Sequences .. " << endl;
    int th = 0;
    uint64_t total = 0, no_kmers = 0;
    int count = 0;
    uint64_t hash=0, fhVal=0, rhVal=0;
    while ((l = kseq_read(seq)) >= 0) {
        ++total;
        //cout << "\r" << total << " processing ..." << flush;
        int len = strlen(seq->seq.s);
        ntHashIterator itr(seq->seq.s, 1, n);
        while (itr != itr.end()) {
            hash = (*itr)[0];
            ++no_kmers;
            uint8_t tz = trailing_zeros(hash);
            if(tz >= th){
                //uint32_t hash = 0; 
                //MurmurHash3_x86_32((uint8_t *)&hash1, 8, 0, &hash); 
                //if(MAP.find(hash) != MAP.end()) MAP[hash].second += 1;  // increment the counter if already there 
                if(MAP[tz].find(hash) != MAP[tz].end()) MAP[tz][hash] += 1; 
                //if(MAP.find(hash) != MAP.end()) MAP[hash].second += 1;
                else{ //// insert if not there 
                  //MAP.insert(make_pair(hash, make_pair(tz, 1))); 
                  MAP[tz].insert(make_pair(hash, 1)); 
                  ++count;  // insert if not there 
                  //cout << "\r" << "count: " << count << flush;// << endl;
                  if(count == k){
                    int cnt = MAP[th].size();
                    count = count - cnt;
                    SMap().swap(MAP[th]);
                    //MAP[th].clear(); //MAP[th].resize(0);
                    ++th;
                    cout  << "count: " << count << endl; 
                   /* int cnt = MAP[th].size();
                    MAP[th].clear();
                    count = count - cnt; 
                    ++th;*/
                   /*for (auto it = MAP.begin(); it != MAP.end(); ++it){
                      if (it->second.first == th) { MAP.erase(it); } //it = MAP.erase(it);
                      //else ++it;
                    }
                    count = MAP.size();
                    ++th;*/
                    //cout << "th: " << th << endl;
                    /*decltype(MAP) newmap;
                    for (auto&& p : MAP)
                        if(p.second.first > th)
                            newmap.emplace(move(p));
                    MAP.swap(newmap);*/
                    //count = MAP.size();
                    //++th;
                  }
                }
            }
	    ++itr;
	}
    }
    //exit(0); 
    cout << "th: " << th << endl;
    cout << "No. of sequences: " << total << endl;
    FILE *fo = fopen(outf.c_str(), "w");
    uint32_t csize = 0; //MAP.size();
    for(int i=th; i<64; i++) csize += MAP[i].size();
    unsigned long F0 = csize * pow(2, (th));
    cout << "F0: " << F0 << endl;
    fprintf(fo, "F1\t%lu\n", no_kmers);
    fprintf(fo, "F0\t%lu\n", F0);
    cout << endl;
    cout << "total: " << total << endl;
    cout << "no_kmer: " << no_kmers << endl;
    //unsigned long freq[65]; 
   unsigned long *freq = new unsigned long[cov];
   for(int i=1; i<=cov; i++) freq[i] = 0;
    unsigned long tot = 0;
    int xx = 0;
    for(int i=th; i<64; i++){
      for(auto& p: MAP[i]){
      //for(auto& p: MAP){
        if(p.second <= cov) freq[p.second]++;
      }
    } 
    /*for (auto it = m.begin(); it != m.end(); it++){
      if(it->second <= 65) freq[it->second]++;
    }*/
    //FILE *fo = fopen(argv[4], "w"); 
    cout << "th: " << th << endl;
    for(int i=1; i<=cov; i++){
      unsigned long fff = (freq[i]*pow(2, th));
      fprintf(fo, "f%d\t%lu\n", i, fff); 
    }
    fclose(fo);
    //unsigned long F0 = MAP.size() * pow(2, (th)); 
    //cout << "F0: " << F0 << endl;
    return 0;
        
}
