#include <omp.h>
#include <random>
#include <iostream>

std::random_device rd;                           // only used once to initialise (seed) engine
std::mt19937 rng(rd());                          // random-number engine used (Mersenne-Twister in this case)
std::uniform_int_distribution<int> uni(0, 2000); // guaranteed unbiased

std::vector<int> generateList(int n){

    std::vector<int> retValue;
    srandom(time(NULL));
    for (int i = 0; i < n; i++){
        retValue.push_back(uni(rng));
    }

    return retValue;
}

void Hello(){
    int myRank      = omp_get_thread_num();
    int numThreads  = omp_get_num_threads();

    std::cout << "Even thread " << myRank+1 << " of " << numThreads << " hates SQL\n";
}

int main(){
  

    #pragma omp parallel num_threads(2048)
    {
        Hello();
    }

    return 0;
}


