// Software: Gaussian Elimination Algorithm 
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
    
    // The matrix dimensions
    int nRows = mxGetM(input_pointers[0]);
    int nCols = mxGetN(input_pointers[0]);
    
    // Get the matrix
    double **A = new double* [nRows];
    for (int i = 0; i < nRows; ++i) {
        A[i] = new double [nCols];
    }
    
    vector2matrix(mxGetPr(input_pointers[0]), nRows, nCols, A);
    
    // Gaussian Elimination
    int determinant_sign = 1;
    int N = min(nRows, nCols);
    bool is_singular = false;
    
    for (int k = 0; k < N; ++k) {
        // Find the k-th pivot
        int i_max = k;
        for (int i = k + 1; i < nRows; ++i) {
            if (abs(A[i][k]) > abs(A[i_max][k])) {
                i_max = i;
            }
        }
        
        // Check if the matrix is singular or not
        if (A[i_max][k] == 0) {
            std::cerr << "The matrix is singular!" << std::endl;
            is_singular = true;
            continue;
        }
        
        // Swap the k-th row and the i_max row
        for (int j = k; j < nCols; ++j) {
            double temp = A[k][j];
            A[k][j] = A[i_max][j];
            A[i_max][j] = temp;
        }
        
        // For compute the determinant
        if (k != i_max) {
            determinant_sign = - determinant_sign;
        }
        
        // Do this for all rows below the pivot
        for (int i = k + 1; i < nRows; ++i) {
            double f = A[i][k] / A[k][k];
            for (int j = k + 1; j < nCols; ++j) {
                A[i][j] -= f * A[k][j];
            }
            A[i][k] = 0.0;
        }
    }
    
    // For computing the determinant
    double determinant = determinant_sign;
    if ((is_singular) || (nRows != nCols)) {
        determinant = 0.0;
    } else {
        for (int i = 0; i < nRows; ++i) {
            determinant *= A[i][i];
        }
    }
    
    // Return the row-echelon form matrix
    output_pointers[0] = mxCreateDoubleMatrix(nRows, nCols, mxREAL);
    matrix2vector(A, nRows, nCols, mxGetPr(output_pointers[0]));
    
    // Return the determinant
    output_pointers[1] = mxCreateDoubleScalar(determinant);
}