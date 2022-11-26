#include <iostream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <fstream>
#include <string.h>
#include <array>

#define Matrix std::vector< std::vector<float> >

// Generates random symmetric matrix
Matrix genRandomMatrix(int n) {
    Matrix A;
    A.resize(n);
    for(int i=0; i<n; ++i) {
        A[i].resize(n);
    };
    for(int i = 0; i<n; ++i) {
        for(int j=i; j<n; ++j) {
            A[i][j] = rand() % 100;
            A[j][i] = A[i][j];
        }
    };
    return(A);
}

void printMatrix(Matrix A) {
    for(int i=0; i<A.size(); ++i) {
        printf("\t");
        for(int j=0; j<A[i].size(); ++j) {
            printf("%.5f ", A[i][j]);
        };
        printf("\n");
    };
}

void saveMatrix(Matrix A, std::string file_name) {
    std::ofstream out(file_name);
    for(int i=0; i<A.size(); ++i) {
        for(int j=0; j<A[i].size(); ++j) {
            out << A[i][j] << " ";
        };
        out << "\n";
    };
    out.close();
}


void JMstep(Matrix &A, int i, int j) {
    float xi = (A[j][j] - A[i][i])/(2*A[i][j]);
    float t = 1./(abs(xi)+sqrt(1.+xi*xi));
    if(xi<0) {
        t = -t;
    };
    float cos = 1./sqrt(1.+t*t);
    float sin = t*cos;
    float oldA_ii = A[i][i];
    float oldA_jj = A[j][j];
    float oldA_ij = A[i][j];
    for(int q=0; q<A.size(); q++) {
        float oldA_iq = A[i][q];
        // float oldA_qi = A[q][i];
        float oldA_jq = A[j][q];
        // float oldA_qj = A[q][j];
        A[i][q] = cos * oldA_iq - sin * oldA_jq;
        A[q][i] = cos * oldA_iq - sin * oldA_jq;

        A[j][q] = sin * oldA_iq + cos * oldA_jq;
        A[q][j] = sin * oldA_iq + cos * oldA_jq;
    };
    A[i][i] = oldA_ii - t*oldA_ij;
    A[j][j] = oldA_jj + t*oldA_ij;
    A[i][j] = 0;
    A[j][i] = 0;
}

double offDiagSum(Matrix A) {
    double sum = 0;
    for(int i=0; i<A.size(); ++i) {
        for(int j=i+1; j<A[i].size(); j++) {
            sum += pow(A[i][j], 2);
        };
    };
    return sum;
}

int getNextI(int k, int n);
int getNextJ(int k, int n);

int getNextI(int k, int n) {
    if(k == 0) return 1;
    int ik = getNextI(k - 1, n);
    int jk = getNextJ(k - 1, n);
    if(ik < n - 1 && jk < n)  return ik;
    if(ik < n - 1 && jk == n) return ik + 1;
    return 1;
}

int getNextJ(int k, int n) {
    if(k == 0) return 2;
    int ik = getNextI(k - 1, n);
    int jk = getNextJ(k - 1, n);
    if(ik < n - 1 && jk < n)  return jk+1;
    if(ik < n - 1 && jk == n) return ik + 2;
    return 2;
}

std::vector< std::array<int, 2> > genIJs(int n) {
    std::vector< std::array<int, 2> > ijs;
    int maxK = n*(n-1)/2;
    ijs.resize(maxK);
    int k=0;
    std::array<int, 2> ij;
    for(int i=0; i<n; i++) {
        for(int j=i+1; j<n; j++) {
            ij[0]=i; ij[1]=j;
            if(k<maxK) {
                ijs[k]=ij;
            }
            else {
                printf("Problem with i=%d, j=%d, k=%d\n", i, j, k);
            }
            k++;
        }
    }
    return ijs;
}

int main(int argc, char **argv) {
    int N = 50;
    if(argc>1) {
        N = atoi(argv[1]);
    };
    Matrix A = genRandomMatrix(N);
    saveMatrix(A, "A.txt");
    double init_sum = offDiagSum(A);
    printf("sum(A) = %f\n", offDiagSum(A));

    std::vector< std::array<int, 2> > ijs = genIJs(N);
    int step = ijs.size()/10;
    std::ofstream progress_file("./progress.txt");
    for(int it = 0; it<3; ++it) {
        printf("== %d ",it);
        for(int k = 0; k<ijs.size(); ++k) {
            if( k % step == 0) {
                std::cout<<"."<<std::flush;
            }
            int i = ijs[k][0];
            int j = ijs[k][1];
            JMstep(A, i, j);
            double offSum = offDiagSum(A);
            progress_file << offSum << "\n";
        };
        std::cout << std::endl;
    }
    progress_file.close();
    double final_sum = offDiagSum(A);
    printf("sum(A)_final  = %f\n", final_sum);
    printf(" ratio = %f\n", final_sum/init_sum);
    saveMatrix(A, "JA_20");


    return 0;
}