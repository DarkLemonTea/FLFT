#include<stdlib.h>
#include<math.h>
#include<string.h>

#include"../FT_CC/tol_multigather.h"
#include"../Ring_FD/ring_FD.h"

#ifndef HEADER_FILE
#define HEADER_FILE
#include"../FT_CC/tol_bcast.h"
#endif

#define Swap(x,y) {double* temp; temp = x; x = y; y = temp;}

/*********************************************************************/
void Gen_matrix(
	double  **A_local  /* out */,
	int    n        /* in 矩阵大小 */,
	int    my_rank  /* in 进程号*/,
	int    comm_size        /* in 进程数 */) {

	int       i, j;
	double     row_sum;
	int       n_bar, rem, new_rem;

	rem = n % comm_size;
	if (rem == 0) {
		n_bar = n / comm_size;
		for (i = 0; i < n_bar; i++) {
			for (j = 0; j < n; j++) {
				if ((my_rank*n_bar + i) == j) {
					A_local[i][j] = (double)(n / 2 + n / 10);
				}
				else
				{
					if ((i % 3 == 0) || (j % 2 == 0)) A_local[i][j] = (double)(0.7);
					else if ((i % 2 == 0) || (j % 3 == 0)) A_local[i][j] = (double)(0.3);
					else A_local[i][j] = (double)(0.5);
				}
			}
		}
	}
	else
	{
		n_bar = n / (comm_size - 1);
		if (my_rank != comm_size - 1) {
			for (i = 0; i < n_bar; i++) {
				for (j = 0; j < n; j++) {
					if ((my_rank*n_bar + i) == j) {
						A_local[i][j] = (double)(n / 2 + n / 10);
					}
					else
					{
						if ((i % 3 == 0) || (j % 2 == 0)) A_local[i][j] = (double)(0.7);
						else if ((i % 2 == 0) || (j % 3 == 0)) A_local[i][j] = (double)(0.3);
						else A_local[i][j] = (double)(0.5);
					}
				}
				for (j = n; j < comm_size*n_bar; j++) {
					A_local[i][j] = 0;
				}
			}
		}
		else {
			new_rem = n % (comm_size - 1);
			for (i = 0; i < new_rem; i++) {
				for (j = 0; j < n; j++) {
					if ((my_rank*n_bar + i) == j) {
						A_local[i][j] = (double)(n / 2 + n / 10);
					}
					else
					{
						if ((i % 3 == 0) || (j % 2 == 0)) A_local[i][j] = (double)(0.7);
						else if ((i % 2 == 0) || (j % 3 == 0)) A_local[i][j] = (double)(0.3);
						else A_local[i][j] = (double)(0.5);
					}
				}
				for (j = n; j < comm_size*n_bar; j++) {
					A_local[i][j] = 0;
				}
			}
			//补零
			for (i = new_rem; i < n_bar; i++) {
				for (j = 0; j < comm_size*n_bar; j++) {
					A_local[i][j] = 0;
				}
			}
		}	
	}
}/*Gen_matrix*/

/*********************************************************************/
void Gen_vector(
	double  *x_local   /* out */,
	int    n          /* in 矩阵大小 */,
	int    my_rank   /* in 进程号*/,
	int    comm_size         /* in 进程数 */) {

	int   i;
	int   n_bar, rem, new_rem;

	rem = n % comm_size;

	if (rem == 0) {
		n_bar = n / comm_size;
		for (i = 0; i < n_bar; i++){
			if (i % 2 == 0) x_local[i] = (double)(3);
			else if (i % 3 == 0) x_local[i] = (double)(5);
			else x_local[i] = (double)(2);
		}		
	}
	else {
		n_bar = n / (comm_size - 1);
		if (my_rank != comm_size - 1) {
			for (i = 0; i < new_rem; i++) {
				if (i % 2 == 0) x_local[i] = (double)(3);
				else if (i % 3 == 0) x_local[i] = (double)(5);
				else x_local[i] = (double)(2);
			}
		}
		else {
			for (i = new_rem; i < n_bar; i++) {
				x_local[i] = 0;
			}
		}
	}
}/*Gen_vector*/
/*********************************************************************/

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

  /*********************************************************************/

double max_dis(double x[], int n) {
	double max_num = 0;
	int i;
	for (i = 0; i < n; i++) {
		if (x[i] > max_num) max_num = x[i];
	}
	return max_num;
}

/*******************************************************************/
void get_send_recv(
	int handle_num,
	void* sendbuf,
	void* databuf,
	int datacount,
	int my_rank,
	int comm_size,
	int size_of_datatype
) {
	int ind;
	char *tmp_buf = NULL, *tmp_data_buf = NULL;

	for (ind = 0; ind < handle_num; ind++) {
		tmp_buf = (char*)sendbuf + ((MPI_Aint)datacount * (MPI_Aint)ind * (MPI_Aint)size_of_datatype);
		tmp_data_buf = (char*)databuf + 
			((MPI_Aint)datacount * (MPI_Aint)((my_rank - handle_num + ind + comm_size) % comm_size) * 
			(MPI_Aint)size_of_datatype);
		memcpy(tmp_buf, tmp_data_buf, datacount * size_of_datatype);
	}
	tmp_buf = (char*)sendbuf + ((MPI_Aint)datacount * (MPI_Aint)handle_num * (MPI_Aint)size_of_datatype);
	tmp_data_buf = (char*)databuf + ((MPI_Aint)datacount * (MPI_Aint)my_rank * (MPI_Aint)size_of_datatype);
	memcpy(tmp_buf, tmp_data_buf, datacount * size_of_datatype);
}

/*******************************************************************/
/*
能整除，均匀分配
不能整除，每个进程多负责一个，最后一个进程负责的数据补零
*/
int tol_parallel_jacobi(
	//intput
	double **A_local,
	double x_local[],
	double b_local[],
	int n,
	double tol,
	int max_iter,
	MPI_Comm comm
) {
	
	int     i_local, i_global, i, j, loc_up_lim;
	int     n_bar, rem;
	int     iter_num;
	double   *x_temp, *send_x, *send_err;
	double   *x_local_old, *x_local_new;
	double   local_dis, global_dis;
	int		res, my_rank, comm_size;
	double	*max_err;
	int handle_num;
	
	int dtysize = sizeof(double);

	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &comm_size);

	rem = n % comm_size;
	if (rem == 0) {
		n_bar = n / comm_size;
	}
	else {
		n_bar = n / (comm_size - 1);
	}

	//初始化故障检测的变量
	Detector sp;
	init_detector(&sp);
	FD_var fd;
	init_fd_var(0.9, 5, 120, 0.2, 0.2, &fd);
	
	/* Compute the local upper limit index: */
	loc_up_lim = n_bar;

	max_err = (double*)calloc(comm_size,sizeof(double));
	x_temp = (double*)calloc(n_bar * comm_size,sizeof(double));
	
	x_local_old = (double*)calloc(n_bar,sizeof(double));
	x_local_new = (double*)calloc(n_bar,sizeof(double));

	//空间分配问题待商榷
	send_x = (double*)calloc(n_bar * comm_size,sizeof(double));
	send_err = (double*)calloc(comm_size,sizeof(double));

	/* Initialize x_temp */
	memcpy(send_x, b_local, n_bar * dtysize);
	tol_bruck_multigather(sp.Lagging_procs, 0, send_x, MPI_DOUBLE, x_temp, MPI_DOUBLE, n_bar, comm);

	//初始化
	memcpy(x_local_old, b_local, n_bar * dtysize);
	memcpy(x_local_new, b_local, n_bar * dtysize);

	/*迭代开始*/
	iter_num = 1;
	while (iter_num < max_iter) {

		//故障模拟
		if ((my_rank == 2) && (iter_num == 50)) sleep(10);
		if ((my_rank == 5) && (iter_num == 10)) sleep(10);

		Swap(x_local_old, x_local_new);

		for (i_local = 0; i_local < loc_up_lim; i_local++) {
			i_global = i_local + my_rank * n_bar;

			x_local[i_local] = b_local[i_local];
			for (j = 0; j < i_global; j++)
				x_local[i_local] = x_local[i_local] - A_local[i_local][j] * x_temp[j];
			for (j = i_global + 1; j < n; j++)
				x_local[i_local] = x_local[i_local] - A_local[i_local][j] * x_temp[j];
			x_local[i_local] = x_local[i_local] / A_local[i_local][i_global];
		}
		//printf("rank %d stage %d local compute finish\n", my_rank, iter_num);

		//for (i = 0; i < n_bar; i++) {
		//	x_temp[n_bar * my_rank + i] = x_local[i];
		//	x_local_new[i] = x_local[i];
		//}
		memcpy(&x_temp[n_bar * my_rank], x_local, n_bar * dtysize);
		memcpy(x_local_new, x_local, n_bar * dtysize);
	
		/*Barrier*/
		res = ring_FD(iter_num, fd, comm, &sp);

		if (res == FD_SUCCESS) {
			//检测是否有从故障中恢复的进程
			probe_revive_procs(iter_num, comm, &sp);
		}
		
		//需要捡回时，才执行该操作
		if (iter_num % 20 == 0 && res != FD_REVIVE && sp.Lagging_procs.num > 0) { 
			ring_retrieve_procs(iter_num, comm, fd.T_retrieve, &sp); 
			//捡回后直接激活
			activate_revive_procs(comm, &sp);
		}
		if (res == FD_REVIVE) {
			//最简单的迭代追赶方式
			printf("rank %d revive!,current stage is %d\n", my_rank, sp.current_stage);
			iter_num = sp.current_stage;
		}

		/*
		multigather
		根据L找出进程需要多负责的数目
		如果handle>0，应先容错填补缺失值，这里直接取旧值
		*/
		handle_num = get_handle_num(sp.Lagging_procs, my_rank, comm_size);
		if (handle_num > 0) {
			printf("rank %d handle_num %d, iter_num %d\n", my_rank, handle_num, iter_num);
		}
		
		get_send_recv(handle_num, send_x, x_temp, n_bar, my_rank, comm_size, dtysize);
		//printf("rank %d send x %f, iter_num %d\n", my_rank, send_x[0], iter_num);

		//容错多聚集
		tol_bruck_multigather(sp.Lagging_procs, handle_num, send_x, MPI_DOUBLE,
			x_temp, MPI_DOUBLE, n_bar, comm);

		local_dis = Distance(x_local_new, x_local_old, n_bar);
		//printf("rank %d local_dis %f, iter_num %d\n", my_rank,local_dis, iter_num);
		max_err[my_rank] = local_dis;
		
		//多聚集误差
		get_send_recv(handle_num, send_err, max_err, 1, my_rank, comm_size, dtysize);
		//printf("rank %d send err %f, iter_num %d\n", my_rank, send_err[0], iter_num);

		tol_bruck_multigather(sp.Lagging_procs, handle_num, send_err, MPI_DOUBLE,
			max_err, MPI_DOUBLE, 1, comm);

		global_dis = max_dis(max_err, comm_size);

		if (my_rank == 0) {
			printf("rank %d global_dis %f,iter_num %d\n\n", my_rank, global_dis, iter_num);
		}
		
		if (global_dis < tol)
			break;

		iter_num++;
	}

	res = ring_FD(iter_num + 1, fd, comm, &sp);
	if (my_rank == 0) {
		printf("global_dis %f,iter_num %d\n", global_dis, iter_num);

		for (i = 0; i < n; i++) printf("%f ", x_temp[i]);
		printf("\n");
	}

	free(send_x);
	free(send_err);
	free(x_local_new);
	free(x_local_old);
	free(x_temp);
	free(max_err);

	if (global_dis < tol) {
		return 1;
	}
	else
		return 0;
} /*Tol Jacobi */
