#include <iostream>
#include <omp.h>
#include "utils.h"

void calc_eigenvalues(int N) {
    int myRank;
    Matrix A;                                   // Matrix to analyze
    std::vector< std::array<int, 2> > ijs;      // Array with (i,j) pairs
    int k;                                      // Loop counter for process number
    int i, j;                                   // indicies of elements to rotate
    int q;                                      // index for horizontal/vertical updates
    float xi, t, sin, cos;                      // rotation parameters 
    float tempAii, tempAjj;
    std::vector<float> tempAiq, tempAqi, tempAjq, tempAqj;
    

    A = genRandomMatrix(N);
    saveMatrix(A, "pA.txt");
    ijs = genIJs(N);
    std::ofstream ijs_out("ijs.txt");
    for(k=0; k<ijs.size(); ++k) {
        ijs_out << k << " " << ijs[k][0] << " " << ijs[k][1] << "\n";
    };
    ijs_out.close();

    printf("[ijs] = %d\n", ijs.size());

    #pragma omp parallel default(none) shared(A, ijs, N) \
    private(myRank, k, i, j, q, \
    xi, t, sin, cos, \
    tempAii, tempAjj, tempAiq, tempAqi, tempAjq, tempAqj)
    {
        myRank = omp_get_thread_num();
        if(myRank==0) {
            printf("Running with %d threads\n", omp_get_num_threads());
        };

        tempAiq.resize(N);
        tempAjq.resize(N);
        tempAqi.resize(N);
        tempAqj.resize(N);

        #pragma omp for
        for(k=0; k<ijs.size(); ++k) {
            i = ijs[k][0];
            j = ijs[k][1];

            // calculating the rotation parameters
            xi = (A[j][j] - A[i][i])/(2*A[i][j]);
            t = 1./(abs(xi)+sqrt(1.+xi*xi));
            if(xi<0) {
                t = -t;
            };
            cos = 1./sqrt(1.+t*t);
            sin = t*cos;
            
            // calculating new elements (storing in temporary vars)
            tempAii = A[i][i] - t*A[i][j];
            tempAjj = A[j][j] + t*A[i][j];
            for(q=0; q<N; ++q) {
                tempAiq[q] = cos * A[i][q] - sin * A[j][q];
                tempAqi[q] = cos * A[i][q] - sin * A[j][q];
                tempAjq[q] = sin * A[i][q] + cos * A[j][q];
                tempAqj[q] = sin * A[i][q] + cos * A[j][q];
            };
            

            // updating the matrix
            #pragma omp critical(matrix_update)
            {
                A[i][j] = 0;
                A[j][i] = 0;
                A[i][i] = tempAii;
                A[j][j] = tempAjj;
                for(q=0; q<N; ++q) {
                    A[i][q] = tempAiq[q];
                    A[j][q] = tempAjq[q];
                    A[q][i] = tempAqi[q];
                    A[q][j] = tempAqj[q];
                }
            }
        };
    } // ending parallel;
    saveMatrix(A, "pA_final.txt");
}

int main(int argc, char **argv) {
    int N = 10;      // number of rows, columns
    int numThreads;  // number of threads

    if(argc>1) {
        N = atoi(argv[1]);
    };
    numThreads = N*(N-1)/2;
    omp_set_num_threads(numThreads);
    calc_eigenvalues(N);    
    return 0;
}