% Software: Gaussian Elimination Algorithm 
% Author: Hy Truong Son
% Position: PhD Student
% Institution: Department of Computer Science, The University of Chicago
% Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
% Website: http://people.inf.elte.hu/hytruongson/
% Copyright 2016 (c) Hy Truong Son. All rights reserved.

function [] = example_solving_linear_system()
    % Size of the matrix
    N = 1024;
	fprintf('N = %d\n', N);
    
    % Randomize a semi-positive definite matrix
    A = rand([N N]);
    A = A' * A;
    
    b = rand([N 1]);
    
    % Solve the linear system by Gauss Elimination
    fprintf('Gauss Elimination\n');
    B = [A, b];
    B = gauss_elimination_cpp(B);
    
    x1 = zeros(N, 1);
    x1(N) = B(N, N + 1) / B(N, N);
    for i = 2 : N
        row = N - i + 1;
        x1(row) = (B(row, N + 1) - B(row, row + 1 : N) * x1(row + 1 : N)) / B(row, row);
    end
        
    % Compute the solution by Matlab
    fprintf('Matlab inverse computation\n');
    x2 = (A^-1) * b;
        
    % Check the result
    fprintf('||A * x1 - b|| = %.6f\n', sum(abs(A * x1 - b)));
    fprintf('||A * x2 - b|| = %.6f\n', sum(abs(A * x2 - b)));
    fprintf('||x1 - x2|| = %.6f\n', sum(abs(x1 - x2)));
end