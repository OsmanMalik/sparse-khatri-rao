function C = sparse_khatrirao(A, B)
% SPARSE_KHATRIRAO Computes the Khatri-Rao product of sparse matrices
%
% C = SPARSE_KHATRIRAO(A) returns the Khatri-Rao product C of the sparse
% matrices stored in the cell A. Note that A must be a row or column cell
% containing sparse matrices, all of which must have the same number of
% columns.
%
% C = SPARSE_KHATRIRAO(A, B) returns the Khatri-Rao product C of the two
% sparse matrices A and B. Note that A and B must have the same number of
% columns.
%
% The latest version of this code is provided at
% https://github.com/OsmanMalik/sparse-khatri-rao 

% Author:   Osman Asif Malik
% Email:    osman.malik@colorado.edu
% Date:     January 5, 2019

%% Check inputs

% One input argument
if nargin == 1
    % Ensure input is cell
    if ~iscell(A) 
        error('When passing a single argument, it must be a cell.');
    end
    
    % Ensure input is 2-dim cell
    if ~ismatrix(A) 
        error('A cell input should be a column or row cell array.');
    end
    
    % Ensure input is vector cell array
    if size(A, 1) > 1 && size(A, 2) 
        error('A cell input should be a column or row cell array.');
    end
    
    % Check sparsity
    N = length(A);
    for n = 1:N
        if ~issparse(A{n})
            error('Input matrices must be sparse');
        end
    end
    
    % Check number of columns
    for n = 1:N-1
        if size(A{n}, 2) ~= size(A{n+1}, 2)
            error('All input matrices must have the same number of columns.');
        end
    end
    
% Two input arguments
elseif nargin == 2
    
    % Check sparsity
    if ~issparse(A) || ~issparse(B)
        error('Input matrices must be sparse');
    end
    
    % Check number of columns
    if size(A, 2) ~= size(B, 2)
        error('All input matrices must have the same number of columns.');
    end
    
% Return error if nargin is not 1 or 2
else
    error('Must pass one or two arguments to spare_khatrirao.');
end

%% Call C code for computation

addpath('help_functions')

if nargin == 1
    C = sparse_khatrirao_c(A);
elseif nargin == 2
    C = sparse_khatrirao_c({A, B});
end

end