#include "utils.h"
#include "math.h"
#include <omp.h>
#include "CStopWatch.h"

using std::cout;
using std::endl;

// makes household transformation of matrix A on column i
void makeHTstep(Matrix A, Matrix &Aout, int col) {
    int N=A.size();
    float nx2 = 0, ny2=0;
    Matrix H, my_H;
    Vector v, my_V;
    float v2=0;

    genZeroMatrix(N, H);
    genZeroMatrix(N, my_H);

    v.resize(N);
    my_V.resize(N);
    genZeroMatrix(N, Aout);

    #pragma omp parallel default(none) \
        shared(A, v, nx2, ny2, col, N, v2, H, Aout) \
        firstprivate(my_V)
    {
        #pragma omp for reduction(+:nx2)
        for(int i=0; i<N; ++i) {
            nx2 += pow(A[i][col], 2);
        };

        #pragma omp for reduction(+:ny2)
        for(int i=0; i<=col; ++i) {
            ny2 += pow(A[i][col], 2);
        };

        #pragma omp single
        {
            v[col+1] = A[col+1][col] - sqrt(nx2-ny2);
        };
        #pragma omp for
        for(int i=col+2; i<N; ++i) {
            v[i] = A[i][col];
        };

        #pragma omp for reduction(+: v2)
        for(int i=0; i<N; ++i) {
            v2 += v[i]*v[i];
        }

        #pragma omp for
        for(int i=0; i<N; ++i) {
            H[i][i] = 1;
        }
        #pragma omp for collapse(2)       
        for(int i=0; i<N; ++i) {
            for(int j=0; j<N; ++j) {
                H[i][j] += -2*v[i]*v[j]/v2;
            };
        };

        #pragma omp for collapse(4)
        for(int i=0; i<N; ++i) {
            for(int j=0; j<N; ++j) {
                for(int i1=0; i1<N; ++i1) {
                    for(int i2=0; i2<N; ++i2) {
                        Aout[i][j] += H[i1][i]*A[i1][i2]*H[i2][j];
                    };
                };
            };
        };
    };
};

// For the 3-diagonal matrix T  of eigenvalues smaller then mu
// see https : //www5 . in . tum . de/lehre/vorlesungen/parnum/WS16/lecture_ 12. pdf
int getNEV(float mu, Matrix T) {
    int nCh=0;
    int n = T.size();

    // printMatrix(T);

    float d = mu - T[0][0];
    for(int j=1; j<n; ++j) {
        if(d>0) {
            nCh++;
        };
        d = (mu - T[j][j])-pow(T[j][j-1], 2)/d;
        // printf("j=%d, dj=%f, nCh=%d\n", j, d, nCh);
    };
    if(d>0) {
        nCh++;
    };
    return nCh;
}


// Calculates i-th eigenvalue of the 3-diagonal matrix T
float calcEV(int i, Matrix T, float a0, float b0, int max_it = 10) {
    float a = a0, b = b0, c = (a+b)/2;
    int wa;
    for(int it=0; it<max_it; ++it) {
        c = (a+b)/2;
        wa = getNEV(c, T);
        // printf("\t a=%f, b=%f, c= %f, nEv=%d\n", a, b, c, wa);
        if(wa >= i) {
            b = c;
        } else {
            a = c;
        };
    };
    return c;
};

int main(int argc, char **argv) {
    // reads para,s
    int N=10;
    int num_th = 10;

    if(argc>1) {
        N = atoi(argv[1]);
    };
    if(argc>2) {
        num_th = atoi(argv[2]);
    }
    omp_set_num_threads(num_th);
    int debug = 0;
    if(argc>3) {
        debug = atoi(argv[3]);
    }


    // Init and save matrix
    Matrix A;                                   // Matrix to analyze
    A = genRandomMatrix(N);
    if(debug) {
        saveMatrix(A, "in_hA.txt");
    };

    CStopWatch timer;
    timer.startTimer();
    // Make the transform
    Matrix Aout;
    for(int col=0; col<N-2; ++col) {
    // for(int col=0; col<1; ++col) {
        makeHTstep(A, Aout, col);
        for(int i=0; i<N; ++i) {
            for(int j=0; j<N; ++j) {
                A[i][j] = Aout[i][j];
            };
        };
    };
    timer.stopTimer();
    printf("%d %d %f\n", N, num_th, timer.getElapsedTime());

    // clear small elements
    for(int i=0; i<N; ++i) {
        for(int j=0; j<N; ++j) {
            if( abs(Aout[i][j])<1e-4) {
                Aout[i][j] = 0;
            }
        }
    };
    if(debug) {
        saveMatrix(Aout, "out_hA.txt");
    };

    float mu = 90;
    int nn = getNEV(mu, Aout);
    cout << " There are " << nn << " eigenvalues smaller then " << mu << endl;
    std::ofstream ev_out("ev_out.txt");
    for(int iev=1; iev<=N; ++iev) {
        float ev = calcEV(iev, Aout, -200, 600);
        ev_out << iev << " " << ev << endl;
    };
    ev_out.close();
    return 0;
}