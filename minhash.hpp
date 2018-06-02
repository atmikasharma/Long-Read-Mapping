#include <iostream>
#include <string>
#include <list>
#include <cmath>
#include <inttypes.h>
#include <algorithm>  
#include <functional>
#include <cctype>
#include <cstring>
#include <cstdlib>
#include <ctime>

//template <typename T >
class CountEstimator
{
	/* A simple bottom n-sketch MinHash implementation.
    n is the number of sketches to keep
    Still don't know what max_prime is */
public:
	CountEstimator();
	CountEstimator(int n, long max_prime, int ksize, std::string input_file_name, 
				char save_kmers, std::list<long> *hash_list, bool rev_comp);
	void parse_file(bool rev_comp);
	void down_sample(long h);
	void add(std::string kmer, bool rev_comp=false);
	void add_sequence(std::string seq, bool rev_comp);
	void jaccard_count();
	float jaccard(CountEstimator& other);
	void common_count();
	long common(CountEstimator& other);
	void _truncate();
	void _export();
	void count_vector();
	void jaccard_vector();
	std::list<std::string>& get_kmers_list(){return _kmers;};
	long countOverlaps(std::list<long>& x1, std::list<long>& x2, long p);

public:
	int n;
	long max_prime = 9999999999971;
	int ksize = 0;
	long p;
	std::string input_file_name = "";
	char save_kmers = 'n';
	std::list<long> *hash_list = NULL;
	bool rev_comp=false;
	std::list<long> _mins ;
	std::list<long> _counts ;
	std::list<std::string> _kmers ;
	long _true_num_kmers;
};


#define tbl \
"                                                                "\
  /*ABCDEFGHIJKLMNOPQRSTUVWXYZ      abcdefghijklmnopqrstuvwxyz    */\
" TVGH FCD  M KN   YSAABW R       TVGH FCD  M KN   YSAABW R"
  //" TVGH FCD  M KA   YSAABWARA      TVGH FCD  M KA   YSAABWARA"

std::string _revcomp(const std::string& kmer);


uint64_t MurmurHash64A ( const void * key, int len, unsigned int seed );
bool is_prime(int number);
long get_prime_lt_x(long t);
