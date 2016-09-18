% Software: Gaussian Jordan Elimination Algorithm (Reduced row-echelon form)
% Author: Hy Truong Son
% Position: PhD Student
% Institution: Department of Computer Science, The University of Chicago
% Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
% Website: http://people.inf.elte.hu/hytruongson/
% Copyright 2016 (c) Hy Truong Son. All rights reserved.

function [A] = gauss_jordan_elimination(A)
    nRows = size(A, 1);
    nCols = size(A, 2);
    N = min(nRows, nCols);
    
    for k = 1 : N
        % Find the k-th pivot
        [~, i_max] = max(abs(A(k : nRows, k)));
        i_max = i_max + k - 1;
        
        if A(i_max, k) == 0
            fprintf('The matrix is singular!');
            return
        end
        
        % Swap the k-th row and the i_max row
        temp = A(k, k : nCols);
        A(k, k : nCols) = A(i_max, k : nCols);
        A(i_max, k : nCols) = temp(:);
        
        % Do for all rows below the pivot
        for i = k + 1 : nRows
            f = A(i, k) / A(k, k);
            % Do for all remaining elements in the current row
            A(i, k + 1 : nCols) = A(i, k + 1 : nCols) - f * A(k, k + 1 : nCols);
            A(i, k) = 0;
        end
    end
    
    % Reduced row-echelon form
    for k = N : -1 : 1
        A(k, k : nCols) = A(k, k : nCols) / A(k, k);
        for i = k - 1 : -1 : 1
            A(i, k : nCols) = A(i, k : nCols) - A(i, k) * A(k, k : nCols);
        end
    end
end