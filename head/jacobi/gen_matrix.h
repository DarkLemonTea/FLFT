#include <stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

int MIN(int a, int b){
	return ((a) < (b)) ? (a) : (b);
}

/*
Ĭ�Ϸ��󣬷���
�������Ⱦ��󣨷����죩T
���ɶԽ�����ֵ����D
���ɽ�����x
����ʽ����A=TDT'
ֵ����m=Ax,������Ľ������x������

�ֿ��������local����
�����M*M
һ��N*N��
�ܹ�ģNM*NM

Diagȫ1
Xȫ2
��������A(i,j)ֵ��Ӧ��pow(x,(i/M)+(j/M))
��k��aȡֵx=pow(-1,k)*(1-5k/(N*N*M))
*/

//���ɶԽ�D
void Gen_diag_block_D(
	int ind,           //�Խ�����
	int N,
	int M,
	double *diag
) {
	int i = 0;
	double m = (double)M;
	for (i = 0; i < M; i++) {
		diag[i] = 1;
	}
}

//
void Gen_matrix_block_A(
	int row,			//������
	int col,			//������
	int N,              //���п���Ŀ
	int M,				//���С(M*M)
	double **matrix		//����
) {
	int scale = N * M / 2;
	int i = 0, j = 0;
	int row_ind, col_ind;
	for (i = 0; i < M; i++) {
		for (j = 0; j < M; j++) {
			row_ind = row * M + i;
			col_ind = col * M + j;
			if ((row_ind <= scale) && (col_ind <= scale)) {
				matrix[i][j] = (double)MIN(scale + row_ind, scale + col_ind);
			}
			else if ((row_ind > scale) && (col_ind > scale)) {
				matrix[i][j] = (double)MIN(3 * scale - row_ind, 3 * scale - col_ind);
			}
			else if((row_ind <= scale) && (col_ind > scale)) {
				matrix[i][j] = (double)(2 * scale + row_ind - col_ind);
			}
			else {
				matrix[i][j] = (double)(2 * scale - row_ind + col_ind);
			}
			if (row_ind != col_ind) {
				matrix[i][j] /= (double)(M);
			}
			else {
				matrix[i][j] *= 0.8 * N;
				if (((row_ind == 0) && (col_ind < scale))||
					((row_ind < scale) && (col_ind == 0))) {
					matrix[i][j] *= 20;
				}
			}
		}
	}
}

void Gen_vector_block_B(
	int row,			//������
	int N,              //���п���Ŀ
	int M,				//���С(M*M)
	double *b			//b����
) {
	int k = 0, i = 0, j = 0;
	int scale = N * M / 2;
	double val = 0;

	for (k = 0; k < M; k++) {
		//��ѭ����ÿһ��b[i],Ϊ��row*M+k�е������
		b[k] = 0;
		i = row * M + k;
		for (j = 0; j < N*M; j++) {
			if ((i <= scale) && (j <= scale)) {
				val = (double)MIN(scale + i, scale + j);
			}
			else if ((i > scale) && (j > scale)) {
				val = (double)MIN(3 * scale - i, 3 * scale - j);
			}
			else if ((i <= scale) && (j > scale)) {
				val = (double)(2 * scale + i - j);
			}
			else {
				val = (double)(2 * scale - i + j);
			}

			if (i != j) {
				b[k] += val / (double)(M);
			}	
			else {			
				if (((i == 0) && (j < scale)) ||
					((i < scale) && (j == 0))) {
					val *= 20;
				}

				b[k] += val * 0.8 * N;
			}
		}	
	}
}

void Gen_iter_matrix_block_T_and_C(
	int row,			//������
	int col,			//������
	int N,              //���п���Ŀ
	int M,				//���С(M*M)
	double **A,
	double *B,
	double **T,		//����
	double *C
) {
	int i = 0, j = 0;
	double aii = 0;
	int diag_ind;
	int scale = N * M / 2;

	for (i = 0; i < M; i++) {
		//������е�aii��Ϊ����
		aii = 0;
		diag_ind = row * M + i;
		if (diag_ind <= scale) {
			aii = scale + diag_ind;
		}
		else {
			aii = 3 * scale - diag_ind;
		}
		if (i == 0) { aii *= 20; }
		aii *= 0.8 * N;

		C[i] = B[i] / aii;

		for (j = 0; j < M; j++) {
			if (row * M + i == col * M + j) {
				T[i][j] = 0;
			}
			else {
				T[i][j] = -A[i][j] / aii;
			}
		}
	}
}

void directly_gen_iter_matrix_block_T_and_C(
	int row,			//������
	int col,			//������
	int N,              //���п���Ŀ
	int M,				//���С(M*M)
	double **T,		//����
	double *C
) {
	int i = 0, j = 0;
	double aii = 0;
	int diag_ind;
	int scale = N * M / 2;

	double **A_local;
	double *b_local;

	A_local = (double**)calloc(M, sizeof(double*));
	for (i = 0; i < M; i++) {
		A_local[i] = (double*)calloc(M, sizeof(double));
	}
	b_local = (double*)calloc(M, sizeof(double));

	Gen_matrix_block_A(row, col, N, M, A_local);

	Gen_vector_block_B(row, N, M, b_local);

	for (i = 0; i < M; i++) {
		//������е�aii��Ϊ����
		aii = 0;
		diag_ind = row * M + i;
		if (diag_ind <= scale) {
			aii = scale + diag_ind;
		}
		else {
			aii = 3 * scale - diag_ind;
		}
		if (i == 0) { aii *= 20; }
		aii *= 0.8 * N;

		C[i] = b_local[i] / aii;

		for (j = 0; j < M; j++) {
			if (row * M + i == col * M + j) {
				T[i][j] = 0;
			}
			else {
				T[i][j] = -A_local[i][j] / aii;
			}
		}
	}

	for (i = 0; i < M; i++) {
		free(A_local[i]);
	}
	free(A_local);
	free(b_local);
}


void Gen_init_value_X(
	int M,
	double value,
	double *x
) {
	int i = 0;
	for (i = 0; i < M; i++) {
		x[i] = value;
	}
}