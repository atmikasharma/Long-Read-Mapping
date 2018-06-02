#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <fstream>
#include "minhash.hpp"
// #include "bloom_filter.hpp"
#include "containmentminhash.cpp"
using namespace std;

bool KmerCompare(string &str1, string &str2){
   return str1.compare(str2);
}

float trueJaccards(vector<string> *aKmerVectB, vector<string> *aKmerVectA) {
   std::vector<string> *v = new vector<string>;
   v->reserve(aKmerVectA->size() + aKmerVectB->size());
   std::vector<string>::iterator it, itEnd;

   sort(aKmerVectB->begin(), aKmerVectB->end());
   sort(aKmerVectA->begin(), aKmerVectA->end());

   // bigger set has to go as first argument
   itEnd=std::set_intersection(aKmerVectB->begin(), aKmerVectB->end(), aKmerVectA->begin(), aKmerVectA->end(), v->begin(), KmerCompare);
   uint countIntersect = 0;
   if (itEnd-v->begin() > 0){
      countIntersect = itEnd-v->begin();
   }
   std::cout << "The intersection has " << countIntersect << " elements:\n";

   itEnd=std::set_union(aKmerVectB->begin(), aKmerVectB->end(), aKmerVectA->begin(), aKmerVectA->end(), v->begin(), KmerCompare);
   uint countUnion = 0;
   if (itEnd-v->begin() > 0){
      countUnion = itEnd-v->begin();
   }
   std::cout << "The union has " << countUnion << " elements:\n";

   float trueJaccard = countIntersect / (float)countUnion;

   delete v;
   return trueJaccard;
}

void generateKmers(string &aLargeString, const uint &aKmerSize, vector<string> *aKmerVect, const uint numKmers){
   vector<string>::iterator it = aKmerVect->begin();
   for (uint i=0;i < numKmers; i++){

      aKmerVect->insert(it, aLargeString.substr(i, aKmerSize));
      it++;
   }
}

void printKmers(vector<string> *aKmerVect){
   vector<string>::iterator it;
   for (it = aKmerVect->begin(); it < aKmerVect->end(); it++){
      cout << *it << ' ';
   }
   cout << endl;
}

int main()
{
   string small_string = "CTACGCAAGGGTACCGCGACCTCAATGCTTGCACGTTACTCCTCTGCCGTTAACATGTCTTCAGGTTGATCAAGCTGACTAGCACTTGATTTCCAGGATA";
   string base_names = "../data/Genomes/base_name.txt";
   // Kmer size
   const uint ksize = 11;
   // probability error rate for bloom filter
   const float prob_error_rate = 0.01;
   // max number of hashes to retain for each string
   const uint h = 10;
   // num kmers formed from each string
   uint size_A = small_string.length() - ksize + 1 ;
    uint size_A_est = 0;
   vector<string> *kmersA = new vector<string>;
   kmersA->reserve(size_A);
   generateKmers(small_string, ksize, kmersA, size_A);
   vector<float> *min_jaccards = new vector<float>;
   CountEstimator ch(10, 9999999999971, 11, "", 'y', NULL, false);
   ch.add_sequence(small_string, false);
   vector<string> *kmersB = new vector<string>;
   string ek;
   while(getline("kmer_list.txt"),line){
     kmersB.push_back(line);
   }
   CountEstimator re();
   string line;
   while(getline("minhash.txt"),line){
     re._mins.insert(line);
   }
   float minestimate = re.jaccard(ch);

   bloom_parameters parameters;
   // How many elements roughly do we expect to insert?
   parameters.projected_element_count = 1.15*large_string.length();
   // Maximum tolerable false positive probability? (0,1)
   parameters.false_positive_probability = prob_error_rate; // 1 in 100
   // Simple randomizer (optional)
   parameters.random_seed = 0xA5A5A5A5;

   parameters.compute_optimal_parameters();
   //Instantiate Bloom Filter
   bloom_filter filter(parameters);

   // Insert kmers into Bloom Filter
   for (size_t i = 0; i < size_B; ++i)
   {
      if (!filter.contains(kmersB->at(i))){
         filter.insert(kmersB->at(i));
         size_B_est++;
      }
      else{
         cout << "string already present in bloom_filter : " << kmersB->at(i) << endl;
      }
   }
   cout << "size_B_est :" << size_B_est << endl;
   float int_est = 0;

   for (size_t i = 0; i < size_A; ++i)
   {
      if (filter.contains(kmersA->at(i)))
      {
         int_est++;
            //cout << "BF contains: " << kmersA->at(i) << endl;
      }
      else{
         cout << "BF doesnot contain" << kmersA->at(i) << endl;
      }
   }

   cout << "int_est : " << int_est << endl;
   int_est -= h * prob_error_rate;
   cout << "int_est (after adjustment): " << int_est << endl;
   float containment_est, jaccard_est, minhash_est, jaccard_est_minhash;
   containment_est = int_est /(float) h;
   jaccard_est = size_A * containment_est / (size_A + size_B_est - size_A * containment_est);
   jaccard_est_minhash = size_A * minestimate / (size_A + size_B_est - size_A * minestimate );
   cout << "Minhash index estimate: " << minestimate << endl;
   cout << "Containment index estimate: " << containment_est << endl;
   cout << "Jaccard index estimate (via the containment approach): " << jaccard_est << endl;
   cout << "Jaccard index estimate (via the Minhash  approach): " << jaccard_est_minhash << endl;


   cout << "MinHash" << minestimate << endl;

   // to wrote the mean of min and containment hashes between the sample and the values in the hdf5 exports.

}
