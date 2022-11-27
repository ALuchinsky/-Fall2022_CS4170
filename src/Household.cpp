#include "utils.h"
#include "math.h"
#include <omp.h>
#include "CStopWatch.h"


// makes household transformation of matrix A on column i
void makeHTstep(Matrix A, Matrix &Aout, int col) {
    int N=A.size();
    float nx2 = 0, ny2=0;
    Matrix H;
    Vector v;
    float v2=0;
    genZeroMatrix(N, H);
    v.resize(N);
    genZeroMatrix(N, Aout);


    for(int i=0; i<N; ++i) {
        nx2 += pow(A[i][col], 2);
        if(i<=col) {
            ny2 += pow(A[i][col], 2);
        };
    };
    for(int i=0; i<N; ++i) {
        if(i<=col) {
            v[i]=0;
        }
        else if(i==col+1) {
            v[i] = A[i][col] - sqrt(nx2-ny2);
        }
        else {
            v[i] = A[i][col];
        }
    };

    for(int i=0; i<N; ++i) {
        v2 += v[i]*v[i];
    }
    for(int i=0; i<N; ++i) {
        for(int j=0; j<N; ++j) {
            if(i==j) {
                H[i][j] += 1;
            } else {
                H[i][j] = 0;
            }
            H[i][j] += -2*v[i]*v[j]/v2;
        };
    }

    for(int i=0; i<N; ++i) {
        for(int j=0; j<N; ++j) {
            Aout[i][j] = 0;
            for(int i1=0; i1<N; ++i1) {
                for(int i2=0; i2<N; ++i2) {
                    Aout[i][j] += H[i1][i]*A[i1][i2]*H[i2][j];
                };
            };
        };
    };
};

int main(int argc, char **argv) {
    // reads para,s
    int N=10;
    int num_th = 10;

    if(argc>1) {
        N = atoi(argv[1]);
    };
    if(argc>2) {
        num_th = atoi(argv[num_th]);
    }

    // Init and save matrix
    Matrix A;                                   // Matrix to analyze
    A = genRandomMatrix(N);
    saveMatrix(A, "in_hA.txt");

    CStopWatch timer;
    timer.startTimer();
    // Make the transform
    Matrix Aout;
    for(int col=0; col<N-2; ++col) {
        makeHTstep(A, Aout, col);
        for(int i=0; i<N; ++i) {
            for(int j=0; j<N; ++j) {
                A[i][j] = Aout[i][j];
            };
        };
    };
    timer.stopTimer();
    printf("N p time\n");
    printf("%d %d %f\n", N, num_th, timer.getElapsedTime());

    // clear small elements
    for(int i=0; i<N; ++i) {
        for(int j=0; j<N; ++j) {
            if( abs(Aout[i][j])<1e-4) {
                Aout[i][j] = 0;
            }
        }
    };
    saveMatrix(Aout, "out_hA.txt");
    return 0;
}