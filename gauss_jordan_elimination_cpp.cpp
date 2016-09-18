// Software: Gauss-Jordan Elimination Algorithm 
// Author: Hy Truong Son
// Position: PhD Student
// Institution: Department of Computer Science, The University of Chicago
// Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
// Website: http://people.inf.elte.hu/hytruongson/
// Copyright 2016 (c) Hy Truong Son. All rights reserved.

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <set>
#include <iterator>
#include <algorithm>
#include <ctime>

#include "mex.h"

using namespace std;

void vector2matrix(double *input, int nRows, int nCols, double **output) {
    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < nCols; ++j) {
            output[i][j] = input[j * nRows + i];
        }
    }
}

void matrix2vector(double **input, int nRows, int nCols, double *output) {
    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < nCols; ++j) {
            output[j * nRows + i] = input[i][j];
        }
    }
}

void mexFunction(int nOutputs, mxArray *output_pointers[], int nInputs, const mxArray *input_pointers[]) {
    if (nInputs != 1) {
        std::cerr << "The number of input parameters must be exactly 1 (the matrix)!" << std::endl;
        return;
    }
    
    // Matrix dimensions
    int nRows = mxGetM(input_pointers[0]);
    int nCols = mxGetN(input_pointers[0]);
    
    // Matrix
    double **A = new double* [nRows];
    for (int i = 0; i < nRows; ++i) {
        A[i] = new double [nCols];
    }
    
    vector2matrix(mxGetPr(input_pointers[0]), nRows, nCols, A);
    
    // Computation
    int N = min(nRows, nCols);
    bool is_singular = false;
    double f, temp;
    int i_max;
    
    for (int k = 0; k < N; ++k) {
        // Find the k-th pivot
        i_max = k;
        for (int i = k + 1; i < nRows; ++i) {
            if (abs(A[i][k]) > abs(A[i_max][k])) {
                i_max = i;
            }
        }
        
        // Check if the matrix is singular or not
        if (A[i_max][k] == 0.0) {
            is_singular = true;
            continue;
        }
        
        // Swap the k-th row with the i_max row
        for (int j = k; j < nCols; ++j) {
            temp = A[k][j];
            A[k][j] = A[i_max][j];
            A[i_max][j] = temp;
        }
        
        // Do this for all rows below the k-th row
        for (int i = k + 1; i < nRows; ++i) {
            f = A[i][k] / A[k][k];
            for (int j = k + 1; j < nCols; ++j) {
                A[i][j] -= f * A[k][j];
            }
            A[i][k] = 0.0;
        }
    }
    
    // Reduced row-echelon form
    for (int k = N - 1; k >= 0; --k) {
        if (A[k][k] == 0.0) {
            continue;
        }
        for (int j = k + 1; j < nCols; ++j) {
            A[k][j] /= A[k][k];
        }
        A[k][k] = 1.0;
        for (int i = 0; i < k; ++i) {
            f = A[i][k];
            for (int j = k; j < nCols; ++j) {
                A[i][j] -= f * A[k][j];
            }
        }
    }
    
    // Return the reduced row-echelon form
    output_pointers[0] = mxCreateDoubleMatrix(nRows, nCols, mxREAL);
    matrix2vector(A, nRows, nCols, mxGetPr(output_pointers[0]));
    
    if (is_singular) {
        std::cerr << "The matrix is singular!" << std::endl;
    }
}