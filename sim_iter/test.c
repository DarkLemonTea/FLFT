#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>

#include"../head/FLFT/FLFT.h"

#include<sys/time.h>
#include<unistd.h>


int main(){
	int res;
	int my_rank, comm_size;

	MPI_Init(NULL, NULL);
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &comm_size);

	//ָ����������ڴ洢�ٵ������б��ָ������б��ָ��
	Detector sp;
	init_detector(TYPE_RING_TREE, comm, 0, 5, &sp);

	//���ϼ��������������÷����ʺ͵ȴ�ʱ������
	FD_var fd;
	init_fd_var(0.90, 5, 5, 0.2, 0.2, &fd);

	struct timeval start, end;
	double cost_time;
	
	//ģ�������̹���
	if (my_rank != 5 && my_rank != 13) {
		FLFT_FD(0, fd, comm, &sp);

		//printf("rank %d lagging num is %d\n", my_rank, sp.Lagging_procs.num);
		
		retrieve_procs(1, comm, 1000, &sp);
		//printf("rank %d finish retrieve\n", my_rank);

		activate_revive_procs(comm, &sp);
		FLFT_FD(1, fd, comm, &sp);
	}
	else {
		sleep(1);
		FLFT_FD(0, fd, comm, &sp);
		FLFT_FD(1, fd, comm, &sp);
	}

	MPI_Finalize();
	return 0;
}

