// # Here we will pre-compute the min hash sketches to speed the simulation. Could parallelize this
// # to make it go faster, but it only needs to be done once, so I don't really care
// import os
// import screed
// import khmer
// import MinHash as MH
// import bz2
// from multiprocessing.dummy import Pool
// import multiprocessing
// from itertools import *
#include <iostream>
#include <fstream>
#include <string>
#include <thread>
#include <list>
#include "linux/limits.h"
#include <experimental/filesystem>
#include <unordered_set>
#include <mutex>

#include "minhash.hpp"

using namespace std;
namespace fs = std::experimental::filesystem;

unsigned int num_cores = thread::hardware_concurrency();
long long prime = 9999999999971; // taking hashes mod this prime
unsigned int ksize = 21;  // k-mer length
unsigned int max_h = 500;  // max number of hashes in sketch
list<string> file_names; 
unsigned int testCount = 0;
mutex listMutex;
mutex minhashMutex;
list<CountEstimator> MHS;

// string filename = "/fatfs/Fall2017Courses/ComputationalBiology/project/MinHashMetagenomics-master/data/Viruses/FileNames.txt";
string relativePath = "/../../../MinHashMetagenomics-master/data/Viruses/";
string filename = "FileNames.txt";


void generateKmersAddtoCE(CountEstimator& MHS, unordered_set<string>& kmers, string& seq, unsigned int aksize, string& agenome){
	for (int i = 0; i < seq.size()-aksize+1; i++){
		string kmer = seq.substr(i, aksize);
		string kmer_rev = _revcomp(kmer);
		// cout << "kmer : " << kmer;
		// cout << "\nkmer_rev: " << kmer_rev;
		if (kmer < kmer_rev){
			kmers.insert(kmer);
			MHS.add(kmer);
			// cout << "\n--------kmer--------\n";	
		}else{
			kmers.insert(kmer_rev);
			MHS.add(kmer_rev);
			// cout << "\n--------kmer_rev--------\n";
		}
		
	}
	MHS._true_num_kmers = kmers.size();
	MHS.input_file_name = fs::path(agenome).filename().string();
	// cout << "kmers.size : " << kmers.size() << "\n";
}

void make_minhash(string& agenome, unsigned int amax_h, unsigned int aprime, unsigned int aksize){
	unordered_set<string> kmers;
	CountEstimator MHS(amax_h, aprime, aksize, "", 'y', NULL, false);
	
	std::ifstream input(agenome);
	if(!input.good()){
		std::cerr << "Error opening '"<<agenome<<"'. Bailing out." << std::endl;
		return;
	}
	std::string line, name, content;
	while( std::getline( input, line ).good() ){
        if( line.empty() || line[0] == '>' ){ // Identifier marker
            if( !name.empty() ){
            	generateKmersAddtoCE(MHS, kmers, content, aksize, agenome);
            	name.clear();
            }
            if( !line.empty() ){
            	name = line.substr(1);
            }
            content.clear();
        } else if( !name.empty() ){
            if( line.find(' ') != std::string::npos ){ // Invalid sequence--no spaces allowed
            	name.clear();
            	content.clear();
            } else {
            	content += line;
            }
        }
    }
    if( !name.empty() ){ // Print out what we read from the last entry
    	// std::cout << "last : " << name << " : " << content.size() << std::endl;
    }

    generateKmersAddtoCE(MHS, kmers, content, aksize, agenome);

    ofstream fos;
//    fos.open("data/Viruses/" + fs::path(agenome).filename().string() + ".Hash21mers.fa" );
    // cout << "open : " << "data/Viruses/"+fs::path(agenome).filename().string() + ".Hash21mers.fa" << endl;
    fos.open("minhash_list.txt");
    for (auto& mins: MHS._mins){
    	fos << ">\n" << mins << "\n";
    }
    fos.close();
    fos.open("kmer_list.txt");
    for (auto& kmers: MHS._kmers){
        fos << ">\n" << kmers << "\n";
    }
    fos.close();

    // cout << "close" << endl;

	// minhashMutex.lock();
	// testCount++ ;
	// minhashMutex.unlock();
} 

// thread will take a file and process it
void threadFunc(list<string>& aFilenames){
	listMutex.lock();
	while(!aFilenames.empty()){
		string curFile = aFilenames.front();
		// cout << this_thread::get_id() << " filename : " << aFilenames.front() << "\n";
		aFilenames.pop_front();
		listMutex.unlock();
		// do all the processing here
		make_minhash(curFile, max_h, prime, ksize);
		listMutex.lock();
	}
	listMutex.unlock();
}

void createThreads(unsigned int aNumThreads, list<string>& aFilenames){
	vector<thread> thrdVect; // twice to avoid wasting time in thread creation
	for (auto i=0; i < aNumThreads*2; i++){
		thrdVect.push_back(thread(threadFunc,std::ref(aFilenames)));
	}
	for (auto& thrd: thrdVect){
		thrd.join();
		// cout << "thread joined : " << thrd.get_id() << "\n";
	}
	cout << "aFilenames.size() : " << aFilenames.size() << "\n";
	// for(auto& f: aFilenames){
	// 	cout << f << "\n";
	// }
}


int main(){
	ifstream fis;
	string line;
	fs::path p = fs::current_path();
	string absPath = p.string() + relativePath + filename; 
	cout << absPath << "\n";
	fis.open(absPath);
	if (fis.is_open()) {
		while ( getline (fis,line) ) {
			// cout << fs::path(line).filename().string() << "\n";
			file_names.push_back(fs::path(line).string());
		}
	}else {
		cout << "Unable to open file\n";
	}
	fis.close();
	// cout << *file_names.begin() << "\n";
	// for_each(file_names.begin(), file_names.end(), [](string &filename){cout << filename << "\n";});
	// cout << "numfiles : " << file_names.size() << "\n";
	// for_each(file_names.begin(), file_names.end(), createThread);
	createThreads(num_cores, file_names);
	cout << "testCount : " << testCount ;
	cout << endl;
	return 0;
}
