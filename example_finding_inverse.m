% Software: Gaussian Jordan Elimination Algorithm 
% Author: Hy Truong Son
% Position: PhD Student
% Institution: Department of Computer Science, The University of Chicago
% Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
% Website: http://people.inf.elte.hu/hytruongson/
% Copyright 2016 (c) Hy Truong Son. All rights reserved.

function [] = example_finding_inverse()
    % Size of the matrix
    N = 512;
	fprintf('N = %d\n', N);
    
    % Randomize a semi-positive definite matrix
    A = rand([N N]);
    A = A' * A;
    
    % Gauss-Jordan Elimination
    fprintf('Gauss-Jordan Elimination\n');
    B = gauss_jordan_elimination_cpp([A, eye(N)]);
    B = B(:, N + 1 : 2 * N);
    
    % Matlab
    fprintf('Matlab\n');
    C = A^-1;
    
    % Comparison
    fprintf('||A * B - I|| = %.6f\n', sum(sum(abs(A * B - eye(N)))));
    fprintf('||B * A - I|| = %.6f\n', sum(sum(abs(B * A - eye(N)))));
    fprintf('||A * C - I|| = %.6f\n', sum(sum(abs(A * C - eye(N)))));
    fprintf('||C * A - I|| = %.6f\n', sum(sum(abs(C * A - eye(N)))));
    fprintf('||B - C|| = %.6f\n', sum(sum(abs(B - C))));
end