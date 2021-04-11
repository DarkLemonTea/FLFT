#include <stdio.h> 
#include <stdlib.h>
#include<mpi.h>
#include"ring_FD.h"

#ifndef HEADER_FILE
#define HEADER_FILE
#include"ring_head.h"
#endif

//随机数
int randNext(int left, int right, MPI_Comm comm)
{
	struct timeval t;
	gettimeofday(&t, NULL);

	int  nameLen, numProcs, myID;
	char processorName[MPI_MAX_PROCESSOR_NAME];
	MPI_Comm_rank(MPI_COMM_WORLD, &myID);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Get_processor_name(processorName, &nameLen);

	srand((unsigned)t.tv_sec + myID * numProcs + nameLen);
	return rand() % (right - left + 1) + left;
}

//模拟工作
void sim_work(int min, int max, int unit, MPI_Comm comm)
{
	unsigned work_time;
	int rand_num = randNext(min, max, comm);
	//printf("rand num is %d\n",rand_num);
	work_time = (unsigned)(rand_num * unit);
	usleep(work_time);
}

//环barrier
int double_ring_barrier(MPI_Comm comm) {
	int my_rank, comm_size;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &comm_size);

	int left, right;
	left = (comm_size + my_rank - 1) % comm_size;
	right = (my_rank + 1) % comm_size;

	if (my_rank > 0) {
		MPI_Recv(NULL, 0, MPI_BYTE, left, 16, comm, MPI_STATUSES_IGNORE);
	}
	MPI_Send(NULL, 0, MPI_BYTE, right, 16, comm);
	if (my_rank == 0) {
		MPI_Recv(NULL, 0, MPI_BYTE, left, 16, comm, MPI_STATUSES_IGNORE);
	}

	if (my_rank > 0) {
		MPI_Recv(NULL, 0, MPI_BYTE, left, 16, comm, MPI_STATUSES_IGNORE);
	}
	MPI_Send(NULL, 0, MPI_BYTE, right, 16, comm);
	if (my_rank == 0) {
		MPI_Recv(NULL, 0, MPI_BYTE, left, 16, comm, MPI_STATUSES_IGNORE);
	}

	return 0;
}

int main()
{
	int my_rank, comm_size;
	MPI_Init(NULL, NULL);
	struct timeval pro_start, pro_end;
	gettimeofday(&pro_start, NULL);

	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &comm_size);

	//初始化故障检测的变量
	Detector sp;
	init_detector(&sp);
	FD_var fd;
	init_fd_var(0.95, 5, 120, 0.2, 0.2, &fd);
	
	int res;
	int iter = 0;

//=============================================================================================================================
//故障模拟迭代
//=============================================================================================================================
	while (iter < 300)
	{
		//故障模拟
		if ((my_rank == 2) && (iter == 30)) sleep(2);
		if ((my_rank == 5) && (iter == 5)) sleep(2);
		
		//sim_work(2, 10, 100, comm);
		//usleep(5 * 1000);

		res = ring_FD(iter, fd, comm, &sp);
		
		if (my_rank == 0) {
			printf("rank %d iter_num %d\n", my_rank, iter);
		}

		if (iter % 2 == 0 && res == FD_SUCCESS) {
			probe_revive_procs(iter, comm, &sp);
		}
		if (iter % 10 == 0 && res == FD_SUCCESS && sp.Lagging_procs.num > 0) {
			ring_retrieve_procs(iter, comm, fd.T_retrieve, &sp);
			activate_revive_procs(comm, &sp);
		}
		if (res == FD_REVIVE) {
			printf("rank %d revive!,current stage is %d\n", my_rank, sp.current_stage);
			iter = sp.current_stage;
		}

		iter++;
	}

//=============================================================================================================================
//无故障迭代
//=============================================================================================================================

	//while (iter < 200)
	//{
	//	//故障模拟
	//	//if ((my_rank == 2) && (iter == 20)) sleep(20);
	//	//if ((my_rank == 5) && (iter == 5)) usleep(2 * 1000);

	//	//sim_work(2, 10, 100, comm);
	//	//usleep(5 * 1000);

	//	//double_ring_barrier(comm);
	//	MPI_Barrier(comm);

	//	if (my_rank == 0) {
	//		printf("rank %d iter_num %d\n", my_rank, iter);
	//	}

	//	iter++;
	//}

//=============================================================================================================================
//统计迭代时间
//=============================================================================================================================
	gettimeofday(&pro_end, NULL);
	//毫秒为单位
	double max_cost;
	double cost_time = (pro_end.tv_sec - pro_start.tv_sec)* 1000.0 +
		(pro_end.tv_usec - pro_start.tv_usec) / 1000.0;

	MPI_Reduce(&cost_time, &max_cost, 1, MPI_DOUBLE, MPI_MIN, 0, comm);

	if (my_rank == 0) {
		printf("run time is %lf\n", max_cost);
	}

	MPI_Finalize();
	return 0;
}

