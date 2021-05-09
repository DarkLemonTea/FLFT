#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#include<string.h>

#include"../head/jacobi/FT_block_jacobi.h"

#include<sys/time.h>
#include<unistd.h>

int main() {
	int my_rank, comm_size;
	double precision = 0.000000001;
	int max_iters = 10000;
	int i, res = 0;

	//初始化
	MPI_Init(NULL, NULL);
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &comm_size);

	double **A_local;
	double *b_local;
	double **T_local;
	double *c_local;
	double *x_local;

	int N = (int)sqrt(comm_size);
	int M = 1000;

	int row = my_rank / N;
	int col = my_rank % N;

	A_local = (double**)calloc(M, sizeof(double*));
	for (i = 0; i < M; i++) {
		A_local[i] = (double*)calloc(M, sizeof(double));
	}
	b_local = (double*)calloc(M, sizeof(double));

	Gen_matrix_block_A(row, col, N, M, A_local);
	//printf("rank %d local A[0][0] is %f\n", my_rank, A_local[0][0]);

	Gen_vector_block_B(row, N, M, b_local);
	//printf("rank %d local b[0] is %f\n", my_rank, b_local[0]);

	T_local = (double**)calloc(M, sizeof(double*));
	for (i = 0; i < M; i++) { 
		T_local[i] = (double*)calloc(M, sizeof(double));
	}
	c_local = (double*)calloc(M, sizeof(double));
	x_local = (double*)calloc(M, sizeof(double));

	Gen_iter_matrix_block_T_and_C(row, col, N, M, A_local, b_local,
		T_local, c_local);

	Gen_init_value_X(M, -1, x_local);

	for (i = 0; i < M; i++) {
		free(A_local[i]);
	}
	free(A_local);
	free(b_local);

	//res = block_jacobi(T_local, c_local, x_local, N, M, precision, max_iters, comm);
	res = FT_block_jacobi(T_local, c_local, x_local, N, M, precision, max_iters, comm);

	/*if (row == 0) {
		printf("rank %d local x is: ", my_rank);
		for (i = 0; i < M; i++) {
			printf("%f ", x_local[i]);
		}
	}*/
	
	/*释放空间*/
	for (i = 0; i < M; i++) {
		free(T_local[i]);
	}
	free(T_local);
	free(c_local);
	free(x_local);

	MPI_Finalize();

	return 0;
}