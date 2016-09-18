% Software: Compare Matlab and C++ performances
% Author: Hy Truong Son
% Position: PhD Student
% Institution: Department of Computer Science, The University of Chicago
% Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
% Website: http://people.inf.elte.hu/hytruongson/
% Copyright 2016 (c) Hy Truong Son. All rights reserved.

function [] = compare_cpp_matlab()
    %% Supporting functions
    function [x] = linear_solution(B)
        n = size(B, 1);
        x = zeros(size(B, 1), 1);
        x(n) = B(n, n + 1) / B(n, n);
        for i = 2 : n
            row = n - i + 1;
            x(row) = (B(row, n + 1) - B(row, row + 1 : n) * x(row + 1 : n)) / B(row, row);
        end
    end
    
    %% Randomization the inputs 
    % Size of the matrix
    N = 1024;
	fprintf('N = %d\n', N);
    
    % Randomize a semi-positive definite matrix
    A = rand([N N]);
    A = A' * A;
    
    b = rand(N, 1);
    
    %% Test 1: Solving system of linear equations
    % Function implemented by Matlab
    fprintf('Function implemented by Matlab\n');
    x_matlab = linear_solution(gauss_elimination([A, b]));
    
    % Function implemented by C++
    fprintf('Function implemented by C++\n');
    x_cpp = linear_solution(gauss_elimination_cpp([A, b]));
    
    % Check the difference
    fprintf('||A * x_matlab - b|| = %.6f\n', sum(abs(A * x_matlab - b)));
    fprintf('||A * x_cpp - b|| = %.6f\n', sum(abs(A * x_cpp - b)));
    
    %% Test 2: Finding the inverse
    % Function implemented by Matlab
    fprintf('Function implemented by Matlab\n');
    A_matlab = gauss_jordan_elimination([A, eye(N)]);
    I_matlab = A_matlab(:, N + 1 : 2 * N);
    
    % Function implemented by C++
    fprintf('Function implemented by C++\n');
    A_cpp = gauss_jordan_elimination_cpp([A, eye(N)]);
    I_cpp = A_cpp(:, N + 1 : 2 * N);
    
    % Check the difference
    fprintf('||A * I_matlab - I|| = %.6f\n', sum(sum(abs(A * I_matlab - eye(N)))));
    fprintf('||A * I_cpp - I|| = %.6f\n', sum(sum(abs(A * I_cpp - eye(N)))));
    fprintf('||I_matlab * A - I|| = %.6f\n', sum(sum(abs(I_matlab * A - eye(N)))));
    fprintf('||I_cpp * A - I|| = %.6f\n', sum(sum(abs(I_cpp * A - eye(N)))));
end