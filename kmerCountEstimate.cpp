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
//#include "dna_test.h"
#include "ntHashIterator.hpp"

using namespace std;
//KSEQ_INIT(gzFile, gzread)
KSEQ_INIT(int, read)
std::map<char, char> mapp = {{'A', 'T'}, {'C', 'G'}, {'G', 'C'}, {'T', 'A'}, {'N', 'N'}};

// Function to reverse string and return reverse string pointer of that
void ReverseConstString(char *str)
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
}

char complement(char n)
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
}   
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

uint32_t murmur3_32(const uint8_t* key, size_t len, uint32_t seed) {
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

unsigned trailing_zeros(uint64_t n) {
    return n ? __builtin_ctzll(n) : -1;
}
 
int main(int argc, char** argv)
{
    
    if(argc == 1){
      cout << argv[0] << " <seq.fa> <kmerLen> <minHeap_Size> <seed>" << endl;
      exit(0);
    } 
    int n = atoi(argv[2]);
    FILE *fp;
    kseq_t *seq;
    int l;
    fp = fopen(argv[1], "r");
    if( fp == Z_NULL){
      cout<<"File: "<<argv[1] << " does not exist" <<endl;
      exit(1);
    }
    seq = kseq_init(fileno(fp));
    
    int k = atoi(argv[3]); // size of maxheap i.e. sample size
    int seed = atoi(argv[4]);

    unordered_map<uint64_t, pair<int, uint64_t>> MAP;
    unordered_map<uint64_t, int> m;
    priority_queue<pair<uint64_t, int>, std::vector<pair<uint64_t, int> >, CompareBySecond> sample;
    cout << "read the Sequences .. " << endl;
    int th = 0;
    
    long long total = 0, no_kmers = 0;
    int count = 0;
    while ((l = kseq_read(seq)) >= 0) {
        total++;
        if(total%1000000 == 0) cout << "\r" << (total) << " completed" << flush;
        string seqs(seq->seq.s);
        int len = seqs.length();
        std::string kmer = seqs.substr(0, n);
        uint64_t hash, fhVal=0, rhVal=0;
        hash = NTPC64(kmer.c_str(), n, fhVal, rhVal);
        for(int i=0; i<(len-n); i++) {
          no_kmers++;
          hash = NTPC64(seqs[i], seqs[i+n], n, fhVal, rhVal);
          int tz = trailing_zeros(hash);
          if(tz >= th){
            if(MAP.find(hash) != MAP.end()) MAP[hash].second += 1;
            else MAP.insert(make_pair(hash, make_pair(1, tz)));
            count ++;
            if(count >= k){
              for(auto& p: MAP)
              if(p.second.second == th)  MAP.erase(p.first);
            }
            th += 1;
          }
          /*if(m.find(hash) != m.end()) m[hash]++;
          else { 
            //int tz = trailing_zeros(hash);
            if(tz >= th){
              pair<uint64_t, int> ret2(hash, tz);
              pair<uint64_t, int> ret1(hash, 1);
              m.insert(ret1);
              sample.push(ret2);
              count++; //if(count == k) exit(0);
              if(count == k){ 
                while(sample.top().second == th){
                  m.erase(sample.top().first);
                  sample.pop();
                  count--;
                }
                th = th+1;
              }
            }
          }*/
        }
    }
    
   /* cout << "total: " << total << endl;
    cout << "no_kmer: " << no_kmers << endl;
    unsigned long freq[65]; for(int i=1; i<=65; i++) freq[i] = 0;
    unsigned long tot = 0;
    int xx = 0;
    for (auto it = m.begin(); it != m.end(); it++){
      if(it->second <= 65) freq[it->second]++;
    }
    cout << "th: " << th << endl;
    for(int i=1; i<=64; i++){
      unsigned long fff = (freq[i]*pow(2, th));
      printf("f%d\t%lu\n", i, fff); 
    }
    unsigned long F0 = m.size() * pow(2, (th)); 
    cout << "F0: " << F0 << endl;*/
    return 0;
        
}
 
