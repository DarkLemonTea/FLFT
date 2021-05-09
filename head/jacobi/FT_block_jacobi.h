#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>

#include "../FLFT/FLFT.h"

#include "gen_matrix.h"
#include "block_jacobi.h"
#include "backup_strategy.h"

#include<sys/time.h>
#include<unistd.h>

#define FAILURE 2
#define BACKUP_X 500
#define BACKUP_ROOT_X 501
#define DATA_RECOVERY 502

/*******************************************************************/

int FT_block_jacobi(
	double **T_local,		//迭代阵T
	double *c_local,		//值向量
	double *x_local,		//解向量
	int N,                  //N*N个进程
	int M,					//每块规模M，矩阵大小(NM)*(NM)	
	double precision,		//迭代精度
	int max_iter,			//最大迭代次数
	MPI_Comm comm
) {
	int     iter_num, i = 0;
	double   *x_local_temp, *x_local_old, *x_local_new;
	double   local_dis = 0, global_dis = 0;
	int		fd_res = 0, reb_res, res, my_rank, comm_size;

	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &comm_size);

	//指针变量，用于存储迟到进程列表、恢复进程列表的指针
	Detector sp;
	//init_detector(TYPE_RING, comm, 0, comm_size, &sp);
	init_detector(TYPE_RING_TREE, comm, 0, 8, &sp);

	//故障检测变量，用于设置放行率和等待时间上限
	FD_var fd;
	init_fd_var(0.95, 5, 30, 2, 2, &fd);

	MPI_Aint seg_count = (MPI_Aint)M;
	int max_reqs = 0;

	//获得备份信息
	Backup my, row_root, col_root;
	back_up_index(my_rank, N, &my);
	back_up_index((my_rank / N) * N + (my_rank / N), N, &row_root);
	back_up_index( (my_rank % N) * N + (my_rank % N), N, &col_root);

	//printf("my_rank: %d, src: %d des: %d\n", my.rank, my.src.rank, my.des.rank);
	//printf("row_root: %d, src: %d des: %d\n", row_root.rank, row_root.src.rank, row_root.des.rank);
	//printf("col_root: %d, src: %d des: %d\n", col_root.rank, col_root.src.rank, col_root.des.rank);


	/* Compute the local upper limit index: */
	int dtsize = sizeof(double);

	x_local_temp = (double*)calloc(M, dtsize);
	x_local_old = (double*)calloc(M, dtsize);
	x_local_new = (double*)calloc(M, dtsize);

	//初始化
	Gen_init_value_X(M, -1, x_local_old);
	memcpy(x_local_new, x_local_old, M * dtsize);

	//对同行进程的备份数据
	double  *x_src_temp, *x_src_old, *x_src_new, *c_src;
	double	**T_src;
	double  src_local_dis, src_global_dis;
	x_src_temp = (double*)calloc(M, dtsize);
	x_src_old = (double*)calloc(M, dtsize);
	x_src_new = (double*)calloc(M, dtsize);
	T_src = (double**)calloc(M, sizeof(double*));
	for (i = 0; i < M; i++) {
		T_src[i] = (double*)calloc(M, sizeof(double));
	}
	c_src = (double*)calloc(M, sizeof(double));
	//备份数据(矩阵直接生成，但是一般情况下需要根据实际情况发送备份的数据)
	directly_gen_iter_matrix_block_T_and_C(my.src.rank/N, my.src.rank%N,
		N, M, T_src, c_src);
	Gen_init_value_X(M, -1, x_src_old);
	memcpy(x_src_new, x_src_old, M * dtsize);

	//对同列root进程的备份数据
	double  *x_root_temp, *x_root_old, *x_root_new, *c_root;
	double	**T_root;	
	double  cr_local_dis, cr_global_dis;
	if (col_root.des.rank == my_rank) {
		x_root_temp = (double*)calloc(M, dtsize);
		x_root_old = (double*)calloc(M, dtsize);
		x_root_new = (double*)calloc(M, dtsize);
		T_root = (double**)calloc(M, sizeof(double*));
		for (i = 0; i < M; i++) {
			T_root[i] = (double*)calloc(M, sizeof(double));
		}
		c_root = (double*)calloc(M, sizeof(double));
		//备份数据(矩阵直接生成，但是一般情况下需要根据实际情况发送备份的数据)
		directly_gen_iter_matrix_block_T_and_C(col_root.rank/N, col_root.rank%N,
			N, M, T_root, c_root);
		Gen_init_value_X(M, -1, x_root_old);
		memcpy(x_root_new, x_root_old, M * dtsize);
	}

	Coll_Tree r_ct, b_ct;//分别用于reduce和bcast的树
	Coll_Tree srcb_ct, crr_ct, crb_ct;//分别用于col_root的reduce和bcast的树

	creat_cube_tree(my_rank, N, 'r', row_root.rank, &r_ct);
	creat_cube_tree(my_rank, N, 'c', col_root.rank, &b_ct);
	
	struct timeval t1, t2;
	double cost_time;
	gettimeofday(&t1, NULL);

	iter_num = 0;
	
	while (iter_num < max_iter) {
		
		if (iter_num == 5 && my_rank == 11) {
			sleep(5);
		}

		Swap(x_local_old, x_local_new);

		//利用迭代矩阵计算出本地的向量T_local*x_local_old=x_local_tmp
		matrix_multi_vector(M, T_local, x_local_old, x_local_temp);

		/*if (iter_num >= 20) {
			printf("rank %d iter num %d: lagging num is %d, last lagging num is %d\n",
				my_rank, iter_num, sp.Lagging_procs.num, sp.Last_lagging_procs.num);
		}*/

		//-----------------------------------------------------------------------------------------
		fd_res = FLFT_FD(iter_num, fd, comm, &sp);
		
		if (fd_res == FD_FAILURE) {
			printf("FD failure\n");
			res = FAILURE;
			break;
		}
		if (fd_res == FD_REVIVE) {
			//被捡回，恢复到和其他进程一致的状态
			printf("rank %d revive!,current stage is %d\n", my_rank, sp.current_stage);
			MPI_Recv(x_local_new, M, MPI_DOUBLE, my.des.rank, 
				DATA_RECOVERY, comm, MPI_STATUSES_IGNORE);
			
			reb_res = rebuild_comm_graph(sp, my_rank, N, 'r', row_root.rank, &r_ct);
			if (reb_res == REBUILD_FAILURE) { res = FAILURE; break; }
			reb_res = rebuild_comm_graph(sp, my_rank, N, 'c', col_root.rank, &b_ct);
			if (reb_res == REBUILD_FAILURE) { res = FAILURE; break; }

			iter_num = sp.current_stage + 1;
			continue;
		}
		if (fd_res == FD_PASS) {
			//建立通信对象
			if (sp.Lagging_procs.num > 0) {
				reb_res = rebuild_comm_graph(sp, my_rank, N, 'r', row_root.rank, &r_ct);
				if (reb_res == REBUILD_FAILURE) { res = FAILURE; break; }
				reb_res = rebuild_comm_graph(sp, my_rank, N, 'c', col_root.rank, &b_ct);
				if (reb_res == REBUILD_FAILURE) { res = FAILURE; break; }
				/*if (iter_num > 0) {
					printf("rank %d, r_ct: p %d, lc %d, rc %d, b_ct: p %d, lc %d, rc %d\n",
						my_rank, r_ct.parent.rank, r_ct.lchild.rank, r_ct.rchild.rank,
						b_ct.parent.rank, b_ct.lchild.rank, b_ct.rchild.rank);
				}*/
			}

			if (my.src.rank != EMPTY) {
				Swap(x_src_old, x_src_new);
				if (in(sp.Lagging_procs, my.src.rank)) {
					reb_res = rebuild_comm_graph(sp, my.src.rank, N, 'c',
						my.src.rank%N + (my.src.rank%N)*N, &srcb_ct);
					if (reb_res == REBUILD_FAILURE) { res = FAILURE; break; }
				}
			}

			if (col_root.des.rank == my_rank) {
				Swap(x_root_old, x_root_new);
				if (in(sp.Lagging_procs, col_root.rank)) {
					reb_res = rebuild_comm_graph(sp, col_root.rank, N, 'r', col_root.rank, &crr_ct);
					if (reb_res == REBUILD_FAILURE) { res = FAILURE; break; }
					reb_res = rebuild_comm_graph(sp, col_root.rank, N, 'c', col_root.rank, &crb_ct);
					if (reb_res == REBUILD_FAILURE) { res = FAILURE; break; }
				}
			}
		}

		//-----------------------------------------------------------------------------------------
		//计算x_local向量
		//如果src在失败进程列表，额外承担reduce任务
		if (my.src.rank != EMPTY) {
			if (in(sp.Lagging_procs, my.src.rank)) {
				matrix_multi_vector(M, T_src, x_src_old, x_src_temp);
				vector_add_vector(M, x_local_temp, x_src_temp);
			}
		}
		
		//自己负责的部分
		/*printf("rank %d, r_ct: p %d, l %d, r %d\n",
			my_rank, r_ct.parent.rank, r_ct.lchild.rank, r_ct.rchild.rank);*/

		res = Reduce_2D(x_local_temp, x_local_new, M, MPI_DOUBLE, MPI_SUM,
			r_ct, comm, seg_count, max_reqs);
		
		//printf("rank %d finish iter %d reduce\n", my_rank, iter_num);

		vector_add_vector(M, x_local_new, c_local);

		/*printf("rank %d, b_ct: root %d, nextsize %d,  p %d, l %d, r %d\n",
			my_rank, b_ct.root, b_ct.nextsize, b_ct.parent.rank, b_ct.lchild.rank, b_ct.rchild.rank);*/

		res = Bcast_2D(x_local_new, M, MPI_DOUBLE, b_ct, comm, seg_count);
			
		//printf("rank %d finish bcast\n", my_rank);

		//col root进程失败，额外承担其任务
		if ((col_root.des.rank == my_rank) &&
			in(sp.Lagging_procs, col_root.rank)) {
			
			res = Reduce_2D(x_root_temp, x_root_new, M, MPI_DOUBLE, MPI_SUM,
				crr_ct, comm, seg_count, max_reqs);
			vector_add_vector(M, x_root_new, c_root);
			res = Bcast_2D(x_root_new, M, MPI_DOUBLE, crb_ct, comm,
				seg_count);
		}

		//备份源进程失败，承担与该进程同列的广播
		if (in(sp.Lagging_procs, my.src.rank)) {
			res = Bcast_2D(x_src_new, M, MPI_DOUBLE, srcb_ct, comm, seg_count);
		}


		//printf("rank %d finish cal x\n", my_rank);
		//-----------------------------------------------------------------------------------------
		//计算距离
		local_dis = Distance(x_local_new, x_local_old, M);
		//printf("rank %d local dis is %f\n", my_rank, local_dis);
		if (in(sp.Lagging_procs, my.src.rank)) {
			src_local_dis = Distance(x_src_new, x_src_old, M);
			local_dis = local_dis > src_local_dis ? local_dis : src_local_dis;
		}

		res = Reduce_2D(&local_dis, &global_dis, 1, MPI_DOUBLE, MPI_MAX,
			r_ct, comm, seg_count, max_reqs);

		res = Bcast_2D(&global_dis, 1, MPI_DOUBLE, b_ct, comm, seg_count);

		//col root进程失败，额外承担其任务
		if ((col_root.des.rank == my_rank) &&
			in(sp.Lagging_procs, col_root.des.rank)) {

			cr_local_dis= Distance(x_root_new, x_root_old, M);
			
			res = Reduce_2D(&cr_local_dis, &cr_global_dis, 1, MPI_DOUBLE, MPI_MAX,
				crr_ct, comm, seg_count, max_reqs);
			
			res = Bcast_2D(&cr_global_dis, 1, MPI_DOUBLE, crb_ct, comm, seg_count);
		}

		//备份源进程失败，承担与该进程同列的广播
		if (in(sp.Lagging_procs, my.src.rank)) {
			res = Bcast_2D(&src_global_dis, 1, MPI_DOUBLE, srcb_ct, comm, seg_count);
		}

		//printf("rank %d finish cal dis\n", my_rank);
		//-----------------------------------------------------------------------------------------	
		if (my_rank == 0) {
			printf("iter num %d global dis is %f\n", iter_num, global_dis);
		}

		if (global_dis<precision || iter_num>max_iter) {
			break;
		}

		/*if (my_rank % 20 == 0 && iter_num >= 3) {
			printf("rank %d finish iter %d finish, global dis is %lf\n",
				my_rank, iter_num, global_dis);
		}*/

		//-----------------------------------------------------------------------------------------
		//备份数据
		if ((!in(sp.Lagging_procs, my.src.rank))
			&& (my.src.rank != EMPTY)) {
			MPI_Irecv(x_src_new, M, MPI_DOUBLE, my.src.rank, BACKUP_X, comm, &my.src.req);
		}
		
		if ((!in(sp.Lagging_procs, col_root.rank))
			&& (col_root.des.rank == my_rank)) {
			//printf("col root des %d is %d\n", col_root.rank, col_root.des.rank);
			MPI_Irecv(x_root_new, M, MPI_DOUBLE, col_root.rank, BACKUP_ROOT_X, comm, &col_root.des.req);
		}
		if (!in(sp.Lagging_procs, my.des.rank)) {
			//printf("rank %d my des is %d\n", my_rank, my.des.rank);
			if (my_rank / N == my_rank % N) {
				MPI_Isend(x_local_new, M, MPI_DOUBLE, my.des.rank, BACKUP_ROOT_X, comm, &my.des.req);
			}
			else {
				MPI_Isend(x_local_new, M, MPI_DOUBLE, my.des.rank, BACKUP_X, comm, &my.des.req);
			}		
		}	
		
		if ((!in(sp.Lagging_procs, my.src.rank))
			&& (my.src.rank != EMPTY)) {
			MPI_Wait(&my.src.req, MPI_STATUSES_IGNORE);
		}
		if ((!in(sp.Lagging_procs, col_root.rank))
			&& (col_root.des.rank == my_rank)) {
			MPI_Wait(&col_root.des.req, MPI_STATUSES_IGNORE);
		}
		if (!in(sp.Lagging_procs, my.des.rank)) {
			MPI_Wait(&my.des.req, MPI_STATUSES_IGNORE);
		}	

		//printf("rank %d finish backup\n", my_rank);
		//-----------------------------------------------------------------------------------------
		//进程捡回
		//需要捡回时，才执行该操作
		if (res == FD_SUCCESS || res == FD_PASS) {
			//检测是否有从故障中恢复的进程
			probe_revive_procs(iter_num, comm, &sp);
		}

		if (iter_num % 20 == 0 && res != FD_REVIVE && sp.Lagging_procs.num > 0) {
			retrieve_procs(iter_num, comm, fd.T_retrieve, &sp);
			//捡回后直接激活
			activate_revive_procs(comm, &sp);
			
			if (in(sp.Revive_procs, my.src.rank)) {
				MPI_Send(x_src_new, M, MPI_DOUBLE, my.src.rank, DATA_RECOVERY, comm);
			}
			if (in(sp.Revive_procs, col_root.rank)
				&& (col_root.des.rank == my_rank)) {
				MPI_Send(x_root_new, M, MPI_DOUBLE, col_root.rank, DATA_RECOVERY, comm);
			}

			//重新调整通信结构
			reb_res = rebuild_comm_graph(sp, my_rank, N, 'r', row_root.rank, &r_ct);
			if (reb_res == REBUILD_FAILURE) { res = FAILURE; break; }
			reb_res = rebuild_comm_graph(sp, my_rank, N, 'c', col_root.rank, &b_ct);
			if (reb_res == REBUILD_FAILURE) { res = FAILURE; break; }
		}

		//-----------------------------------------------------------------------------------------
		iter_num += 1;	
	}

	//这里可以free一下

	gettimeofday(&t2, NULL);
	cost_time = (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;

	if (my_rank == 0) {
		printf("cost time is %lf\n", cost_time);
	}

	memcpy(x_local, x_local_new, M * dtsize);

	if (res == FAILURE) {
		return res;
	}

	if (global_dis < precision) {
		return ITERATION_CONVERGENCE;
	}
	else {
		return ITERATION_NOT_CONVERGENCE;
	}
	
}
