#include <iostream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <fstream>
#include <string.h>
#include <array>

#define Matrix std::vector< std::vector<float> >
#define Vector std::vector<float>

void genZeroMatrix(int n, Matrix &M);

void printVector(std::vector<float> v, std::string name);

// Generates random symmetric matrix
Matrix genRandomMatrix(int n);

// Prints matrix A to the standard output
void printMatrix(Matrix A);

// saves matrix A to a file
void saveMatrix(Matrix A, std::string file_name);

// Calculates the sum of squared matrix elements on the upper-right triangle
double offDiagSum(Matrix A);


// Returns an array of (i,j) pairs for a matrix size n
std::vector< std::array<int, 2> > genIJs(int n);

