#include "utils.h"

int main(int argc, char **argv) {
    int N=10;
    Matrix A;                                   // Matrix to analyze
    std::vector<float> qIn, qOut;
    double sum;
    float lambda;

    if(argc>1) {
        N = atoi(argv[1]);
    };
    A = genRandomMatrix(N);
    saveMatrix(A, "poA_in.txt");
    qIn.resize(N);
    std::fill(qIn.begin(), qIn.end(), 0);
    qIn[0]=1;

    qOut.resize(N);

    for(int it=0; it<5; ++it) {
        for(int i=0; i<N; ++i) {
            sum = 0;
            #pragma omp parallel for reduction(+: sum)
            for(int j=0; j<N; j++) {
                sum += A[i][j]*qIn[j];
            };
            qOut[i] = sum;
        }
        lambda = qOut[3]/qIn[3]*N;
        sum = 0;
        #pragma omp parallel for reduction(+: sum)
        for(int j=0; j<N; j++) {
            sum += qOut[j]*qOut[j];
        };
        sum = sqrt(sum);
        #pragma omp parallel for
        for(int i=0; i<N; ++i) {
            qOut[i] = qOut[i]/sum;
            qIn[i] = qOut[i];
        };
        lambda = lambda/N;
    };
    std::cout << "qOut=";
    for(int i=0; i<N; ++i) {
        std::cout << qOut[i] << " ";
    };
    std::cout << " lambda = " << lambda << std::endl;
    std::cout << "\n";
    return 0;
}