#include "utils.h"
#include "math.h"

// generates household matrix I - v vT/(vT v)
void genHouse(Vector v, Matrix &P) {
    int N=v.size();
}

// makes household transformation of matrix A on column i
void makeHTstep(Matrix A, Matrix &Aout, int col) {
    int N=A.size();
    float nx2 = 0, ny2=0;
    Matrix H;
    for(int i=0; i<N; ++i) {
        nx2 += pow(A[i][col], 2);
        if(i<=col) {
            ny2 += pow(A[i][col], 2);
        };
    };
    // printf(" nx2=%f, ny2=%f\n", nx2, ny2);
    Vector v;
    v.resize(N);
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
    // printVector(v, "v = ");
    genZeroMatrix(N, H);

    float v2=0;
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


    // genHouse(v, H);
    // printMatrix(H);

    genZeroMatrix(N, Aout);
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
    int N=10;
    Matrix A;                                   // Matrix to analyze
    if(argc>1) {
        N = atoi(argv[1]);
    };
    A = genRandomMatrix(N);
    saveMatrix(A, "in_hA.txt");
    Matrix Aout;

    for(int col=0; col<N-2; ++col) {
        makeHTstep(A, Aout, col);
        for(int i=0; i<N; ++i) {
            for(int j=0; j<N; ++j) {
                A[i][j] = Aout[i][j];
            };
        };
    };
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