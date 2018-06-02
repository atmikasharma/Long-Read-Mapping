#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <sys/types.h>
#include <dirent.h>

void read_directory(const std::string& name, stringvec& v)
{
    DIR* dirp = opendir(name.c_str());
    struct dirent * dp;
    while ((dp = readdir(dirp)) != NULL) {
        v.push_back(dp->d_name);
    }
    closedir(dirp);
}

void sim(int num_genomes, int reads){
  std::vector<std::string> files;

std::vector<int> vecOfRandomNums(num_genomes - 1);
std::generate(vecOfRandomNums.begin(), vecOfRandomNums.end(), []() {
  return rand() % 100;
});

std::vector<std::string> names;
for (unsigned int i = 0; i < files.size(),i++){

}
}


int main(){
  unsigned int num_cores = std::thread::hardware_concurrency();
  int num_genomes = 20;
  int num_reads = 10000;
  int num_replicates = 16;
  double prime = 9999999999971;
  float false_pos = 0.001;
  int ksize = 11;
  int max_hashes = 5000;
}
