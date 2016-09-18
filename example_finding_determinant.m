% Software: Gaussian Elimination Algorithm 
% Author: Hy Truong Son
% Position: PhD Student
% Institution: Department of Computer Science, The University of Chicago
% Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
% Website: http://people.inf.elte.hu/hytruongson/
% Copyright 2016 (c) Hy Truong Son. All rights reserved.

function [] = example_finding_determinant()
    % Size of the matrix
    N = 64;
	fprintf('N = %d\n', N);
    
    % Randomize a semi-positive definite matrix
    A = rand([N N]);
    % A = A' * A;
    
    % Gauss Elimination
    fprintf('- Gauss Elimination:\n');
    
    [~, det1] = gauss_elimination_cpp(A);
    
    fprintf('Determinant: %.6f\n', det1);
    
    % Matlab
    fprintf('- Matlab\n');
    
    det2 = det(A);
    
    fprintf('Determinant: %.6f\n', det2);
    
    % Comparison
    fprintf('|det1 - det2| = %.6f\n', abs(det1 - det2));
end