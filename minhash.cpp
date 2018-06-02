#include "minhash.hpp"


using namespace std;

CountEstimator::CountEstimator()
{
	cout<<"CountEstimator::CountEstimator - Default"<<endl;
}

CountEstimator::CountEstimator(int n, long max_prime, int ksize, string input_file_name,
	char save_kmers, list<long> *hash_list, bool rev_comp)
{
	// cout<<"CountEstimator::CountEstimator - Begin"<<endl;
	if(n == 0)
	{
		cout<<"n=0 exception\n";
		exit(0);
	}
	if(ksize == 0)
	{
		cout<<"ksize=0 exception\n";
		exit(0);
	}
	if((ksize % 2) == 0)
	{
		cout<<"Due to an issue with khmer, only odd ksizes are allowed\n";
		exit(0);
	}
	this->ksize = ksize;
	this->hash_list = hash_list;
	// get a prime to use for hashing
	long p = get_prime_lt_x(max_prime);
	this->p = p;
	// initialize sketch to size n
	this->_mins.assign(n,p);
	// initialize the corresponding counts
	this->_counts.assign(n,0);
	// initialize the list of kmers used, if appropriate

	if(save_kmers == 'y') 	{
		this->_kmers.assign(n,"");
		this->save_kmers = save_kmers;
	}
	else 	{
		this->_kmers.resize(0);
	}
	// Initialize file name (if appropriate)
	this->input_file_name = input_file_name;
	if(this->input_file_name != "\0")
	{
		this->parse_file(rev_comp=rev_comp);
	}
	// Optional container for the true number of k-mers in the genome used to populate the sketch
	this->_true_num_kmers = 0;
	//seed for random number generator
	srand(time(0));
	// cout<<"CountEstimator::CountEstimator - end"<<endl;
}

//template <typename T>
void CountEstimator::parse_file(bool rev_comp)
{
	/* opens a file and populates the CountEstimator with it */

}

//template <typename T>
void CountEstimator::down_sample(long h)
{
	/* This will down-sample a sketch to have exactly h elements
        :param h: number of elements you wish to save
        :return: None */
	this->_mins.resize(h);
	this->_counts.resize(h);
	this->_kmers.resize(h);
}

void CountEstimator::add(string kmer, bool rev_comp)
{
	// cout<<"CountEstimator::add - begin"<<endl;
	/* Add kmer into sketch, keeping sketch sorted, update counts accordingly */
	_mins = this->_mins;
	_counts = this->_counts;
	_kmers = this->_kmers;
	uint64_t h;
	if(rev_comp)
	{
		uint64_t h1 = MurmurHash64A(kmer.c_str(),ksize,rand()); // Have to check the 2nd and 3rd arguments of murmurHash function call
		uint64_t h2 = MurmurHash64A(kmer.c_str(),ksize,rand()); // have to implement "h2 = khmer.hash_murmur3(khmer.reverse_complement(kmer))"
		h=(h1<h2)?h1:h2;
		if(h == h2)
		{
			// kmer = khmer.reverse_complement(kmer)
		}
	}
	else
	{
		h = MurmurHash64A(kmer.c_str(),ksize,rand());
		// cout << "h : " << h << endl;
	}

	h = h % this->p;
	//cout << "h % this->p :" << h << endl;
	if(this->hash_list != NULL && !this->hash_list->empty())
	{
		bool flag = false;
		for(list<long>::iterator it =this->hash_list->begin(); it!=this->hash_list->end(); it++)
		{
			if(h==*it)
			{
				flag = true;
				break;
			}
		}
		if(flag==false)
		{
			cout << "flag==false" << endl;
			return;
		}
	}
	if(h>=*(--_mins.end()))
	{
		// cout << "h>=*(_mins.end())" << endl;
		return;
	}
	// bisect function
	long i=0;
	bool foundIndex = false;
	auto it = lower_bound(this->_mins.begin(), this->_mins.end(), h);
	// if (it != this->_mins.end()) {
	// 	foundIndex = true;
	// }

	i = distance (this->_mins.begin(), it);
	//cout << "i:" << i << "\n" ;

	if (*it == h ){
		auto itrCounts = this->_counts.begin();
		advance(itrCounts, i);
		//cout << "*itrCounts : " << *itrCounts << "\n" ;
		*itrCounts += 1;
		return;
	}
	else {
		auto itrMins = this->_mins.begin();
		advance(itrMins, i);
		_mins.insert( itrMins, h);
		_mins.pop_back();
		// cout << "mins list\n";
		// for (auto v : this->_mins)
		// 	std::cout << v << " ";
		// cout << "\n";
		auto itrCounts = this->_counts.begin();
		advance(itrCounts, i);
		_counts.insert( itrCounts, 1);
		_counts.pop_back();
		// cout << "counts list\n";
		// for (auto v : this->_counts)
		// 	std::cout << v << " ";
		// cout << "\n";
		// if(!_kmers.empty())
		if (this->save_kmers == 'y')
		{
			auto itrkmers = this->_kmers.begin();
			advance(itrkmers, i);
			_kmers.insert(itrkmers, kmer);
			_kmers.pop_back();
			// cout << "kmers list\n";
			// for (auto v : this->_kmers)
			// 	std::cout << v << " ";
			// cout << "\n";
			//cout << "_kmers.size() : " << _kmers.size() << endl;
		}
		return;
	}
	// cout<<"CountEstimator::add - end"<<endl;
}

//template <typename T>
void CountEstimator::add_sequence(string seq, bool rev_comp)
{
	cout<<"CountEstimator::add_sequence - Begin"<<endl;
	/* Sanitize and add a sequence to the sketch. */
	transform(seq.begin(), seq.end(), seq.begin(), ptr_fun<int, int>(toupper));
	for (string::size_type l = 0; l < seq.length(); ++l)
	{
		if(seq[l] != 'A' && seq[l] != 'C' && seq[l] != 'T' && seq[l] != 'G')
		{
			seq[l] = 'G';
		}
	}
	uint numKmers = seq.length() - this->ksize + 1;
	for (uint i=0;i < numKmers; i++){
		add(seq.substr(i, this->ksize), rev_comp);
	}
//	for(list<string>::iterator it =this->_kmers.begin(); it!=this->_kmers.end(); it++)
//        {
//		cout<<*it<<endl;
//	}
	cout<<"CountEstimator::add_sequence - end"<<endl;
}

float CountEstimator::jaccard(CountEstimator& other){
	long truelen = this->_mins.size();
	auto itr = this->_mins.end();
	while (truelen > 0 &&  *--itr == this->p){truelen--;};
	// cout << "truelen : " << truelen << endl;
	if (truelen == 0){
		throw "jaccard ValueError";
	}
	return (float)this->common(other)/truelen;
}

//template <typename T>
void CountEstimator::jaccard_count()
{

}

//template <typename T>
void CountEstimator::common_count()
{

}

long CountEstimator::countOverlaps(list<long>& x1, list<long>& x2, long p){
	auto i = x1.begin(), j = x2.begin();
	long processed = 0, common = 0;
	while ( i != x1.end() && j != x2.end()) {
		// cout << "processed : " << processed << endl;
		while (*i < *j){
			i++;
			if (i == x1.end())
				break;
		}
		while (*i > *j){
			j++;
			if (j == x2.end())
				break;
		}
		if (*i == *j){
			if (*i != p){
				common++;
			}
			i++;
			j++;
		}
		processed++;
	}
	// cout << "common : " << common << endl;
	return common;
}
//template <typename T>
long CountEstimator::common(CountEstimator& other){
	// Calculate number of common k-mers between two sketches.
	if (this->ksize != other.ksize){
		throw "different k-mer sizes - cannot compare";
	}
	if (this->p != other.p){
		throw "different primes - cannot compare";
	}
	return countOverlaps(this->_mins, other._mins, this->p);
}

//template <typename T>
void CountEstimator::_truncate()
{

}

//template <typename T>
void CountEstimator::_export()
{

}

//template <typename T>
void CountEstimator::count_vector()
{

}

//template <typename T>
void CountEstimator::jaccard_vector()
{

}



bool is_prime(int number)
{
	if(number < 2)
	{
		return false;
	}
	if(number == 2)
	{
		return true;
	}
	if((number % 2) == 0)
	{
		return false;
	}
	for(int j=3; j<(int(pow(number,0.5))+1); j+=2)
	{
		if((number%j) == 0)
		{
			return false;
		}
	}
	return true;
}

long get_prime_lt_x(long t)
{
	if(t == 1)
	{
		return 1;
	}
	int i = int(t);
	if( (i%2) == 0 )
	{
		i -= 1;
	}
	while(i > 0)
	{
		if(is_prime(i))
		{
			return i;
		}
		i -=2;
	}
	if(i <= 0)
	{
		cout<<"unable to find a prime number < "<< t <<endl;
		exit(0);
	}
}

//typedef unsigned __int64 uint64_t;

// 64-bit hash for 64-bit platforms
uint64_t MurmurHash64A ( const void * key, int len, unsigned int seed )
{
	const uint64_t m = 0xc6a4a7935bd1e995;
	const int r = 47;
	uint64_t h = seed ^ (len * m);
	const uint64_t * data = (const uint64_t *)key;
	const uint64_t * end = data + (len/8);
	while(data != end)
	{
		uint64_t k = *data++;
		k *= m;
		k ^= k >> r;
		k *= m;
		h ^= k;
		h *= m;
	}
	const unsigned char * data2 = (const unsigned char*)data;
	switch(len & 7)
	{
		case 7: h ^= uint64_t(data2[6]) << 48;
		case 6: h ^= uint64_t(data2[5]) << 40;
		case 5: h ^= uint64_t(data2[4]) << 32;
		case 4: h ^= uint64_t(data2[3]) << 24;
		case 3: h ^= uint64_t(data2[2]) << 16;
		case 2: h ^= uint64_t(data2[1]) << 8;
		case 1: h ^= uint64_t(data2[0]);
		h *= m;
	};
	h ^= h >> r;
	h *= m;
	h ^= h >> r;
	return h;
}

std::string _revcomp(const std::string& kmer)
{
	std::string out = kmer;

	auto from = out.begin();
	auto to = out.end();

	char c;
	for (to--; from <= to; from++, to--) {
		c = tbl[(int)*from];
		*from = tbl[(int)*to];
		*to = c;
	}
	return out;
}


#ifdef _TEST_

void test_jaccard_1(){
	CountEstimator E1(1, 31, 21, "", 'n', NULL, false);
	CountEstimator E2(1, 31, 21, "", 'n', NULL, false);
	list<long> l1 = {1, 2, 3, 4, 5}, l2 = {1, 2, 3, 4, 6};
	E1._mins = l1;
	E2._mins = l2;
	cout << "test_jaccard_1 : " << (E1.jaccard(E2) == (float)4/5) ? true : false;
	cout << endl;
	cout << "test_jaccard_1 : " << (E2.jaccard(E1) == (float)4/5) ? true : false;
	cout << endl;
}

void test_jaccard_2_difflen(){
	CountEstimator E1(1, 31, 21, "", 'n', NULL, false);
	CountEstimator E2(1, 31, 21, "", 'n', NULL, false);
	list<long> l1 = {1, 2, 3, 4, 5}, l2 = {1, 2, 3, 4};
	E1._mins = l1;
	E2._mins = l2;
	cout << "test_jaccard_1 : " << (E1.jaccard(E2) == (float)4/5) ? true : false;
	cout << endl;
	cout << "test_jaccard_1 : " << (E2.jaccard(E1) == (float)4/4) ? true : false;
	cout << endl;
}

void test_yield_overlaps(){
	list<long> x1 = {1, 3, 5}, x2 = {2, 4, 6};
	// cout << "x1 size : " << x1.size();
	cout << "test_yield_overlaps : " << (countOverlaps(x1, x2, 9) == 0) ? true : false ;
	cout << endl;
}

void test_yield_overlaps_2(){
	list<long> x1 = {1, 3, 5}, x2 = {1, 2, 4, 6};
	cout << "test_yield_overlaps_2 : " << (countOverlaps(x1, x2, 9) == 1) ? true : false ;
	cout << endl;
	cout << "test_yield_overlaps_2 : " << (countOverlaps(x2, x1, 9) == 1) ? true : false ;
	cout << endl;
}
void test_yield_overlaps_3(){
	list<long> x1 = {1, 3, 6}, x2 = {1, 2, 6};
	cout << "test_yield_overlaps_3 : " << (countOverlaps(x1, x2, 9) == 2) ? true : false ;
	cout << endl;
	cout << "test_yield_overlaps_3 : " << (countOverlaps(x2, x1, 9) == 2) ? true : false ;
	cout << endl;
}

int main()
{
	string p = "CTTGACCTCGTCATTTCACTTTCACTGCGGAGTTTCGGGCAACCGATGGAAA";
    //uint64_t h = MurmurHash64A(p.c_str(),11,11);
    //cout<<"murmur : "<<h<<endl;
    //list<long> hash_list;
	// CountEstimator ch(5000, 9999999999971, 11, "", 'y', NULL, false);
	// ch.add_sequence(p, false);
	test_jaccard_1();
	test_yield_overlaps();
	test_yield_overlaps_2();
	test_yield_overlaps_3();
	return 0;
}
#endif
