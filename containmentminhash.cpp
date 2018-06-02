#include <iostream>
#include <string>
#include "bloom_filter.hpp"

using namespace std;

// 
class ContainmentMinHash{
public:
	ContainmentMinHash(unsigned long long int aElementCount, double aFalsePostProb, unsigned long long int aRndSeed);
	void add(string& aKmer);
	void addSequence(string& aSeq);
	inline unsigned long long int element_count() const;
	template <typename T>
	inline bool contains(const T& t) const;

private:
	ContainmentMinHash();
	bloom_filter filter;
	// bloom_parameters parameters;
};

ContainmentMinHash::ContainmentMinHash(unsigned long long int aElementCount, double aFalsePostProb, unsigned long long int aRndSeed){
	bloom_parameters parameters;
	parameters.projected_element_count = 1.15*aElementCount;
	parameters.false_positive_probability = aFalsePostProb; // 1 in 100
	parameters.random_seed = aRndSeed;
	parameters.compute_optimal_parameters();
	bloom_filter tmpFilter(parameters);
	this->filter = tmpFilter;
}

void ContainmentMinHash::add(string& aKmer) {
	if (!this->filter.contains(aKmer)) {
		this->filter.insert(aKmer);
	}
	else{
		cout << "string already present in bloom_filter : " << aKmer << "\n";
	}
}

void ContainmentMinHash::addSequence(string& aSeq){
// need to populate if its required
}

inline unsigned long long int ContainmentMinHash::element_count() const {
	return this->filter.element_count();
}

template <typename T>
inline bool ContainmentMinHash::contains(const T& t) const{
	return this->filter.contains(t);
}

#ifdef _TEST_
int main() {
	string p = "CTTGACCTCGTCATTTCACTTTCACTGCGGAGTTTCGGGCAACCGATGGAAA";
	ContainmentMinHash cmh(1,2,3);
	return 0;
}
#endif