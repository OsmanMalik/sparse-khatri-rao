# Sparse Khatri-Rao product
An efficiently implemented Matlab function for computing the Khatri-Rao product of sparse matrices.

## Installing and running
All that needs to be done is compiling the C code, which can be done by navigating to the **help_functions** folder in Matlab and then running **mex sparse_khatrirao_c.c**. The function is then used by calling the m file **sparse_khatrirao** in the root folder.

## Example
The function **sparse_khatrirao** can be called either with two sparse matrices as input, or with a single input in the form of a row or column cell containing sparse matrices. Both cases are illustrated below.

```matlab
A = cell(1,4);
for n = 1:4
	A{n} = sprand(10, 10, .5);
end

% Provide two matrices as input
C1 = sparse_khatrirao(A{1},A{2});

% Provide a cell of matrices as input
C2 = sparse_khatrirao(A);
```

## Author contact information
Please feel free to contact me at any time if you have any questions or would like to provide feedback on this code. I can be reached at osman.malik@colorado.edu.