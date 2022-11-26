#include "utils.h"

#ifndef _UTILS_H
#define _UTILS_H

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

void printVector(std::vector<float> v, std::string name="") {
    std::cout << name << " ";
    for(int i=0; i<v.size(); ++i) {
        printf("%f ", v[i]);
    };
    printf("\n");
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

double offDiagSum(Matrix A) {
    double sum = 0;
    for(int i=0; i<A.size(); ++i) {
        for(int j=i+1; j<A[i].size(); j++) {
            sum += pow(A[i][j], 2);
        };
    };
    return sum;
}

void genZeroMatrix(int n, Matrix &M) {
    M.resize(n);
    for(int i=0; i<n; ++i) {
        M[i].resize(n);
        for(int j=0; j<n; ++j) {
            M[i][j]=0;
        }
    };
}


#endif