#include <iostream>
#include <string>
#include <algorithm>  
#include <vector>

#include "bloom_filter.hpp"
using namespace std;

void generateKmers(string &aLargeString, const unsigned int &aKmerSize, vector<string> *aKmerVect, const unsigned int numKmers){
   vector<string>::iterator it = aKmerVect->begin();
   for (int i=0;i < numKmers; i++){
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

bool KmerCompare(string &str1, string &str2){
   return str1.compare(str2);
}

// check A in B
unsigned int trueJaccards(vector<string> *aKmerVectB, vector<string> *aKmerVectA) {
   std::vector<string> *v = new vector<string>;
   v->reserve(aKmerVectA->size() + aKmerVectB->size());
   std::vector<string>::iterator it, itEnd;

   sort(aKmerVectB->begin(), aKmerVectB->end());
   sort(aKmerVectA->begin(), aKmerVectA->end());

   // bigger set has to go as first argument
   itEnd=std::set_intersection(aKmerVectB->begin(), aKmerVectB->end(), aKmerVectA->begin(), aKmerVectA->end(), v->begin(), KmerCompare);
   unsigned int countIntersect = 0;
   if (itEnd-v->begin() > 0){
      countIntersect = itEnd-v->begin();
   }
   std::cout << "The intersection has " << countIntersect << " elements:\n";

   itEnd=std::set_union(aKmerVectB->begin(), aKmerVectB->end(), aKmerVectA->begin(), aKmerVectA->end(), v->begin(), KmerCompare);
   unsigned int countUnion = 0;
   if (itEnd-v->begin() > 0){
      countUnion = itEnd-v->begin();
   }
   std::cout << "The union has " << countUnion << " elements:\n";

   float trueJaccard = countIntersect / (float)countUnion;
   cout << "trueJaccard : " << trueJaccard << endl;
   delete v;
   return trueJaccard;
}


int main()
{
   //"ABCDEFGHIJKLMNOPQRSTUVWXYZ";
   string small_string = "CTACGCAAGGGTACCGCGACCTCAATGCTTGCACGTTACTCCTCTGCCGTTAACATGTCTTCAGGTTGATCAAGCTGACTAGCACTTGATTTCCAGGATA";
   string large_string = "CCGCATCGACAAGCAGGATCTGGATCTATTTCTCTCTTAAATCCATGTAAGGGACGGCAGAAACCTGCTCCTTCTACTTGCTACATCTTCTAGGGTAGAACGAGACCAGAGCCGTTACTGCGATATGAAATCAGTACCGAACGTTGGAACTTATTCAGTTTTAACCCGGTCCCCGTCGCCCAAATCGGGCTATATCATACCCCCGGGCCAAGTGTACAAGTGCATCGATTAAATGCACTAACGGCGAAAGTAAATGATGGACTTTCCAAGCCTGAGGTGGTAAACGCACTTGAATAGAGTCGACAAATTATCGGCTGACGATGCCTTGTAGACCAGCTTTAACACATGACCAGTATAGACGAGGCGGAACTAAGCAATCCCAAGTTTTCGTGCGAGCTGAAGGACCCGGCTCCACGAGATAGAGCTTGTGTTAACAAGAGGCCTCCGGCTGGAAAGATTGGTGGAAACGGCTGCTGTCACGTTTGCATCTTACCGGATGTGCCCCAATGAGGAGTTGATGAACTGGCTGTGACGCAATGGCGAAGAGGAAACGTCTGTATGGCGGATGTAACGTTTTTGCAACACTCCTCCACAACTGCTCCTTTAAGATGACCATCACGAAAATGAAGCTCGTTCGAAATCTTCAAAGATCCGGGGTATAATTGCGCTTCCGGGAGAAGGCCATATGCGATAGCGGTAAGTTTCCACAGCGTATCCAAAAGCGGAGCTTTACGATCTCCCCAGTAAACTGGCTTGTGTCAAGCGGCGAACCCGAATTTCGACGAACCTAGATATTCTCTGGCGACTAACTACTATGCGGATGGGCCTATTCGGGGGATTCAGCCCGCGATACTAGAGCGTAATTAGCCTCGCAAGAATCTAGGTAGCCCCAAAATAGCTTGCTAAAGCGCTAGGGTGCACTGCAGGCAAAATCGAGGTGACTGTACCCCGAGCCATGCATATAACTGGGGGGTACCCTTCCAATAATTGTTATCATACCATCTGCATAGACATATTTAACGGCTCAGTAAATTCGTCGCCATGCGACCTCCAGCATGATCGGTGGCACTCCGTTGTGCGCGGACTGTGTAAACCGCACG" + small_string;
   
   const unsigned int ksize = 11;

   unsigned int size_B = large_string.length() - ksize + 1 ;
   unsigned int size_A = small_string.length() - ksize + 1 ;

   cout << "size_B :" << size_B << endl;
   cout << "size_A :" << size_A << endl;

   vector<string> *kmersB = new vector<string>;
   kmersB->reserve(size_B);
   generateKmers(large_string, ksize, kmersB, size_B);
   //printKmers(kmersB);
   
   vector<string> *kmersA = new vector<string>;
   kmersA->reserve(size_A);
   generateKmers(small_string, ksize, kmersA, size_A);
   //printKmers(kmersA);


   float trueJaccard = trueJaccards(kmersB, kmersA);

   bloom_parameters parameters;
   // How many elements roughly do we expect to insert?
   parameters.projected_element_count = 1090;
   // Maximum tolerable false positive probability? (0,1)
   parameters.false_positive_probability = 0.01; // 1 in 100
   // Simple randomizer (optional)
   parameters.random_seed = 0xA5A5A5A5;

   parameters.compute_optimal_parameters();
   //Instantiate Bloom Filter
   bloom_filter filter(parameters);

   // Insert into Bloom Filter
   {
      // Insert some strings
      for (size_t i = 0; i < size_B; ++i)
      {
         filter.insert(kmersB->at(i));
      }

   }

   // Query Bloom Filter
   {
      // Query the existence of strings
      for (size_t i = 0; i < size_A; ++i)
      {
         if (filter.contains(kmersA->at(i)))
         {
            //cout << "BF contains: " << kmersA->at(i) << endl;
         }
         else{
            cout << "BF doesnot contain" << kmersA->at(i) << endl;
         }
      }

      string invalid_str_list[] = { "AbCX", "iJkX", "XYZX" };

      // Query the existence of invalid strings
      for (size_t i = 0; i < (sizeof(invalid_str_list) / sizeof(string)); ++i)
      {
         if (filter.contains(invalid_str_list[i]))
         {
            cout << "BF falsely contains: " << invalid_str_list[i] << endl;
         }
      }

      // Query the existence of invalid numbers
      for (int i = -1; i > -100; --i)
      {
         if (filter.contains(i))
         {
            cout << "BF falsely contains: " << i << endl;
         }
      }
   }
   delete kmersB;
   delete kmersA;
   return 0;
}
