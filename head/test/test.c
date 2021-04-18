#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>
#include"ring_FD.h"
#include<sys/time.h>
#include<unistd.h>
//#include"tol_bcast.h"

int main()
{
	int res;
	int my_rank, comm_size;

	MPI_Init(NULL, NULL);
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &comm_size);

	//ָ����������ڴ洢�ٵ������б��ָ������б��ָ��
	Detector sp;
	init_detector(&sp);

	//���ϼ��������������÷����ʺ͵ȴ�ʱ������
	FD_var fd;
	init_fd_var(0.95, 5, 120, 0.2, 0.2, &fd);

	struct timeval start, end;
	double cost_time;
	
	////ģ��5�Ž��̹���
	//if (my_rank != 5) {
	//	ring_FD(0, fd, comm, &sp);
	//	printf("rank %d lagging num is %d\n", my_rank, sp.Lagging_procs.num);
	//	ring_retrieve_procs(1, comm, 1000, &sp);
	//	activate_revive_procs(comm, &sp);

	//	ring_FD(1, fd, comm, &sp);
	//}
	//else {
	//	sleep(1);
	//	ring_FD(0, fd, comm, &sp);
	//	ring_FD(1, fd, comm, &sp);
	//}

	////ģ�������̹���
	//if (my_rank != 5 && my_rank != 8 && my_rank != 13) {
	//	ring_FD(0, fd, comm, &sp);
	//	printf("rank %d lagging num is %d\n", my_rank, sp.Lagging_procs.num);
	//	ring_retrieve_procs(1, comm, 1000, &sp);
	//	activate_revive_procs(comm, &sp);

	//	ring_FD(1, fd, comm, &sp);
	//}
	//else {
	//	sleep(1);
	//	ring_FD(0, fd, comm, &sp);
	//	ring_FD(1, fd, comm, &sp);
	//}

	//ģ�������̹���
	if (my_rank != 5 && my_rank != 13) {
		if (my_rank == 8) {
			usleep(200000);
		}
		ring_FD(0, fd, comm, &sp);
		printf("rank %d lagging num is %d\n", my_rank, sp.Lagging_procs.num);
		ring_retrieve_procs(1, comm, 1000, &sp);
		activate_revive_procs(comm, &sp);

		ring_FD(1, fd, comm, &sp);
	}
	else {
		sleep(1);
		ring_FD(0, fd, comm, &sp);
		ring_FD(1, fd, comm, &sp);
	}

	MPI_Finalize();
	return 0;
}

