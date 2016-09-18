% Software: Gaussian Elimination Algorithm 
% Author: Hy Truong Son
% Position: PhD Student
% Institution: Department of Computer Science, The University of Chicago
% Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
% Website: http://people.inf.elte.hu/hytruongson/
% Copyright 2016 (c) Hy Truong Son. All rights reserved.

function [A, d] = gauss_elimination(A)
    nRows = size(A, 1);
    nCols = size(A, 2);
    N = min(nRows, nCols);
    
    d = 1.0;
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
        
        % For computing the determinant
        if k ~= i_max
            d = - d;
        end
            
        % Do for all rows below the pivot
        for i = k + 1 : nRows
            f = A(i, k) / A(k, k);
            % Do for all remaining elements in the current row
            A(i, k + 1 : nCols) = A(i, k + 1 : nCols) - f * A(k, k + 1 : nCols);
            A(i, k) = 0;
        end
    end
    
    % Compute the determinant
    d = prod(diag(A)) / d;
    
    % Reduced row-echelon form
%     for k = N : -1 : 1
%         A(k, k : nCols) = A(k, k : nCols) / A(k, k);
%         for i = k - 1 : -1 : 1
%             A(i, k : nCols) = A(i, k : nCols) - A(i, k) * A(k, k : nCols);
%         end
%     end
end