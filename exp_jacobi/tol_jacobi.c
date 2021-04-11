#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>
#include"ring_FD.h"
#include"tol_bcast.h"
#include"jacobi.h"

int main() {
	int my_rank, comm_size;
	double tol = 0.01;
	int max_iter = 10000;
	int i;

	//初始化
	MPI_Init(NULL, NULL);
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &comm_size);

	double **A_local;
	double *b_local;
	double *x_local;

	int n_bar = 20;

	A_local = (double**)calloc(n_bar,sizeof(double*));
	for (i = 0; i < n_bar; i++) { A_local[i] = (double*)calloc(n_bar*comm_size, sizeof(double)); }
	b_local = (double*)calloc(n_bar,sizeof(double));
	x_local = (double*)calloc(n_bar,sizeof(double));

	Gen_matrix(A_local, n_bar * comm_size, my_rank, comm_size);
	Gen_vector(b_local, n_bar * comm_size, my_rank, comm_size);

	tol_parallel_jacobi(A_local, x_local, b_local, n_bar * comm_size, tol, max_iter, comm);

	/*释放空间*/
	for (i = 0; i < n_bar; i++) {
		free(A_local[i]);
	}
	free(A_local);
	free(b_local);
	free(x_local);

	MPI_Finalize();

	return 0;
}