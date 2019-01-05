/* 
 * SPARSE_KHATRIRAO_C.C
 *
 * Compute the Khatri-Rao product of a cell containing sparse matrices.
 *
 * C = sparse_khatrirao_c(A) returns the Khatri-Rao product C of the sparse 
 * matrices stored in the cell A.
 *
 * The latest version of this code is provided at
 * https://github.com/OsmanMalik/sparse-khatri-rao
 *
 * There are no safety checks in this C code. Consider using the Matlab
 * wrapper function provided in the link above.
 *
 * Please compile by running "mex sparse_khatrirao_c.c" in Matlab.
 *
 * */
 
 /*
  * Author:	Osman Asif Malik
  * Email: 	osman.malik@colorado.edu
  * Date: 	January 5, 2019
  *
  * */

#include <stdio.h>
#include "mex.h"

/* Declare global variables */
double **a, *b;
mwIndex **a_ir, **a_jc, *b_ir, *b_jc, *a_no_rows, b_no_rows, no_cols, cnt;
mwSize N;

/* Define function which recursively computes column in output matrix */
void compute_output_column(mwIndex c, mwIndex n, double x, mwIndex ind) {
	double x_new;
	mwIndex i, ind_new;
	for(i = a_jc[n][c]; i < a_jc[n][c+1]; ++i) {
		x_new = x*a[n][i];
		ind_new = ind*a_no_rows[n] + a_ir[n][i];
		if(n < N-1) {
			compute_output_column(c, n+1, x_new, ind_new);
		} else {
			b[cnt] = x_new;
			b_ir[cnt] = ind_new;
			++cnt;
		}
	}
}

/* mex interface */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, 
	const mxArray *prhs[]) {
	/* Declare other variables */
	mwSize c, n, b_nnz;

	/* Get input variables */
	N = mxGetDimensions(prhs[0])[1];
	a = malloc(N*sizeof(double *));
	a_ir = malloc(N*sizeof(mwIndex *));
	a_jc = malloc(N*sizeof(mwIndex *));
	a_no_rows = malloc(N*sizeof(mwIndex));
	for(n = 0; n < N; ++n) {
		a[n] = mxGetPr(mxGetCell(prhs[0], n));
		a_ir[n] = mxGetIr(mxGetCell(prhs[0], n));
		a_jc[n] = mxGetJc(mxGetCell(prhs[0], n));
		a_no_rows[n] = mxGetM(mxGetCell(prhs[0], n));		
	}
	no_cols = mxGetN(mxGetCell(prhs[0], 1));
		
	/* Compute no rows in output matrix */
	b_no_rows = 1;
	for(n = 0; n < N; ++n) {
		b_no_rows *= a_no_rows[n];
	}
	
	/* Compute nnz in output matrix */
	b_nnz = 1;
	for(c = 0; c < no_cols; ++c) {
		mwIndex prod = 1;
		for(n = 0; n < N; ++n){
			prod *= a_jc[n][c+1] - a_jc[n][c];
		}
		b_nnz += prod;
	}
	
	/* Create sparse output matrix */
	plhs[0] = mxCreateSparse(b_no_rows, no_cols, b_nnz, mxREAL);
	b = mxGetPr(plhs[0]);
	b_ir = mxGetIr(plhs[0]);
	b_jc = mxGetJc(plhs[0]);
	
	/* Compute jc for output matrix */
	b_jc[0] = 0;
	for(c = 0; c < no_cols; ++c) {
		mwIndex prod = 1;
		for(n = 0; n < N; ++n) {
			prod *= a_jc[n][c+1] - a_jc[n][c];
		}
		b_jc[c+1] = b_jc[c] + prod;
	}
	
	/* Compute non-zero elements and ir vector for output matrix */
	cnt = 0;
	for(c = 0; c < no_cols; ++c) {
		compute_output_column(c, 0, 1.0, 0);
	}
	
	/* Free dynamically allocated memory */
	free(a_no_rows);
	free(a_jc);
	free(a_ir);
	free(a);
}