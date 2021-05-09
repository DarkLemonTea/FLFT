#include<stdlib.h>
#include<math.h>
#include<string.h>

#include<sys/time.h>
#include<unistd.h>

#ifndef HEADER_FILE
#define HEADER_FILE
#include "gen_matrix.h"
#endif

#include "block_parallel.h"

#define ITERATION_CONVERGENCE 0
#define ITERATION_NOT_CONVERGENCE 1

#define Swap(x,y) {double* temp; temp = x; x = y; y = temp;}

double Distance(double x[], double y[], int n) {
	int i;
	double max = 0;
	double diff_abs;

	for (i = 0; i < n; i++) {
		diff_abs = ((x[i] - y[i] > 0) ? x[i] - y[i] : y[i] - x[i]);
		if (diff_abs > max) max = diff_abs;
	}
	return max;
} /* Distance */

//��������double�ͷ���
void matrix_multi_vector(int M, double **T, double *x, double *res_x) {
	int i = 0, j = 0;
	for (i = 0; i < M; i++) {
		res_x[i] = 0;
		for (j = 0; j < M; j++) {
			res_x[i] += T[i][j] * x[j];
		}
	}
}

//������a�ӵ�x��
void vector_add_vector(int M, double *x, double *a) {
	int i = 0;
	for (i = 0; i < M; i++) {
		x[i] = x[i] + a[i];
	}
}

/*******************************************************************/
int block_reduce_to_diag_then_bcast_vector(
	double *x_local_temp,
	double *x_local_new,
	double *C_local,
	int N,
	int M,
	int comm
) {
	//����N*N���ͷֲ�����reduce����bcast
	int	res, my_rank, comm_size;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &comm_size);

	MPI_Aint seg_count = (MPI_Aint)M;
	int max_reqs = 0;
	
	//�����root����
	int row_root = (my_rank / N) * N + (my_rank / N);
	int col_root = (my_rank % N) * N + (my_rank % N);

	Coll_Tree r_ct, b_ct;//�ֱ�����reduce��bcast����
	creat_cube_tree(my_rank, N, 'r', row_root, &r_ct);
	creat_cube_tree(my_rank, N, 'c', col_root, &b_ct);

	res = Reduce_2D(x_local_temp, x_local_new, M, MPI_DOUBLE, MPI_SUM,
		r_ct, comm, seg_count, max_reqs);

	//printf("V rank %d Reduce res is %d\n", my_rank, res);

	vector_add_vector(M, x_local_new, C_local);

	res = Bcast_2D(x_local_new, M, MPI_DOUBLE, b_ct, comm, seg_count);

	//printf("V rank %d Bcast res is %d\n", my_rank, res);

	return res;
}

int block_reduce_to_diag_then_bcast_dis(
	double local_dis,
	double *global_dis,
	int N,
	int M,
	int comm
) {
	//����N*N���ͷֲ�����reduce����bcast
	int	res, my_rank, comm_size;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &comm_size);

	MPI_Aint seg_count = (MPI_Aint)M / 2;
	int max_reqs = 0;

	//�����root����
	int row_root = (my_rank / N) * N + (my_rank / N);
	int col_root = (my_rank % N) * N + (my_rank % N);

	Coll_Tree r_ct, b_ct;//�ֱ�����reduce��bcast����
	creat_cube_tree(my_rank, N, 'r', row_root, &r_ct);
	creat_cube_tree(my_rank, N, 'c', col_root, &b_ct);

	res = Reduce_2D(&local_dis, global_dis, 1, MPI_DOUBLE, MPI_MAX,
		r_ct, comm, seg_count, max_reqs);

	//printf("D rank %d Reduce global dis is %f\n", my_rank, *global_dis);

	res = Bcast_2D(global_dis, 1, MPI_DOUBLE, b_ct, comm, seg_count);

	//printf("D rank %d Bcast global dis is %f\n", my_rank, *global_dis);

	return res;
}

/*******************************************************************/
int block_jacobi(
	double **T_local,		//������T
	double *c_local,		//ֵ����
	double *x_local,		//������
	int N,                  //N*N������
	int M,					//ÿ���ģM�������С(NM)*(NM)	
	double precision,		//��������
	int max_iter,			//����������
	MPI_Comm comm
) {
	int     iter_num;
	double   *x_local_temp, *x_local_old, *x_local_new;
	double   local_dis = 0, global_dis = 0;
	int		res, my_rank, comm_size;

	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &comm_size);

	/* Compute the local upper limit index: */
	int dtsize = sizeof(double);

	x_local_temp = (double*)calloc(M, dtsize);

	x_local_old = (double*)calloc(M, dtsize);
	x_local_new = (double*)calloc(M, dtsize);

	//��ʼ��
	Gen_init_value_X(M, -1, x_local_old);
	memcpy(x_local_new, x_local_old, M * dtsize);

	/*
	�ֿ鲢�е���
	��ÿ������ӵ��һ���������T_local
	��ÿ�����̼����������x_local
	�����õ��������������ص�����T_local*x_local_old=x_local_tmp
	��ͬ�еĽ���,�ԶԽǽ���Ϊroot��reduce(x_local_tmp,op=sum)���õ�x_local_new,
	��ͬ�н��̣��ԶԽǽ���Ϊroot��bcast(x_local_new)��ʹ�ø����̵õ����ε�����������x_local_new
	
	���������local_err=dis(x_local_new,T_local)
	��ͬ�н���allreduce����õ�global_err

	��global<precision||iter_num>max_iters�˳�����
	*/
	struct timeval t1, t2;
	double cost_time;
	gettimeofday(&t1, NULL);

	iter_num = 0;
	while (iter_num < max_iter) {
		Swap(x_local_old, x_local_new);

		//���õ��������������ص�����T_local*x_local_old=x_local_tmp
		matrix_multi_vector(M, T_local, x_local_old, x_local_temp);

		res = block_reduce_to_diag_then_bcast_vector(x_local_temp, x_local_new, 
			c_local, N, M, comm);

		local_dis = Distance(x_local_new, x_local_old, M);
		//printf("rank %d local dis is %f\n", my_rank, local_dis);

		res = block_reduce_to_diag_then_bcast_dis(local_dis, &global_dis, N, M, comm);
		
		//printf("rank %d global dis is %f\n", my_rank, global_dis);

		if (my_rank == 0) {
			printf("iter num %d global dis is %f\n", iter_num, global_dis);
		}

		if (global_dis<precision || iter_num>max_iter) {
			break;
		}	

		iter_num += 1;

	}

	gettimeofday(&t2, NULL);
	cost_time = (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;

	if (my_rank == 0) {
		printf("cost time is %lf\n", cost_time);
	}

	memcpy(x_local, x_local_new, M * dtsize);

	if (global_dis < precision) {
		return ITERATION_CONVERGENCE;
	}
	else {
		return ITERATION_NOT_CONVERGENCE;
	}
}