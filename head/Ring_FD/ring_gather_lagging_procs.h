#ifndef HEADER_FILE
#define HEADER_FILE
#include "ring_head.h"
#endif

//���ھۼ�lagging����
int multigather_lagging_procs(
	MPI_Comm comm,
	Ring ring,
	//output
	int *ring_lagging_num,
	int **ring_lagging_procs
) {
	int my_rank, comm_size;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &comm_size);

	int ide, i;
	if (ring.left_proc.rank > my_rank) { ide = 0; }//��ʼ����
	else if (my_rank > ring.right_proc.rank) { ide = 2; }//ĩβ����
	else { ide = 1; }

	int recv_num = 0, local_num = 0, send_num = 0, total_num = 0;
	int *recv_procs = NULL, *local_procs = NULL, *send_procs = NULL, *total_procs = NULL;
	
	statistics_lagging_procs(ring, my_rank, comm_size, &local_num, &local_procs);

	Comm_proc left_proc, right_proc, ver_proc;
	init_proc(ring.left_proc.rank, &left_proc);
	init_proc(ring.right_proc.rank, &right_proc);

	switch (ide)
	{
	case 0:
	{	//��һ�ֻ�״�ಥ��������ʵ��	
		while (1) {
			send_procs_template(comm, &right_proc, RING_GATHER_LAGGING_NUM, RING_GATHER_LAGGING_PROCS, local_num, local_procs);
			if (right_proc.comm_stage == FINISH) break;
		}

		//�ڶ��ְ��ջ���
		left_proc.comm_stage = READY; right_proc.comm_stage = READY;
		while (1)
		{
			//��������
			if (left_proc.comm_stage != FINISH) {
				recv_procs_template(comm, &left_proc, RING_MULTICAST_LAGGING_NUM, RING_MULTICAST_LAGGING_PROCS, &total_num, &total_procs);
				if (left_proc.comm_stage == FINISH) {
					*ring_lagging_num = total_num;
					*ring_lagging_procs = total_procs;
				}
			}
			else {
				//��������
				send_procs_template(comm, &right_proc, RING_MULTICAST_LAGGING_NUM, RING_MULTICAST_LAGGING_PROCS, total_num, total_procs);
			}
			if (right_proc.comm_stage == FINISH) break;
		}
		break;
	}
	case 1:
	{
		//��һ��	
		while (1) {
			if (left_proc.comm_stage != FINISH) {
				recv_procs_template(comm, &left_proc, RING_GATHER_LAGGING_NUM, RING_GATHER_LAGGING_PROCS, &recv_num, &recv_procs);
			}
			else {
				if (right_proc.comm_stage == READY) {
					merge_procs(local_num, recv_num, local_procs, recv_procs, &send_num, &send_procs);
					send_procs_template(comm, &right_proc, RING_GATHER_LAGGING_NUM, RING_GATHER_LAGGING_PROCS, send_num, send_procs);
				}
				else {
					send_procs_template(comm, &right_proc, RING_GATHER_LAGGING_NUM, RING_GATHER_LAGGING_PROCS, send_num, send_procs);
				}
			}
			if (right_proc.comm_stage == FINISH) break;
		}

		//�ڶ��ְ��ջ���
		left_proc.comm_stage = READY; right_proc.comm_stage = READY;
		while (1)
		{
			//��������
			if (left_proc.comm_stage != FINISH) {
				recv_procs_template(comm, &left_proc, RING_MULTICAST_LAGGING_NUM, RING_MULTICAST_LAGGING_PROCS, &total_num, &total_procs);
				if (left_proc.comm_stage == FINISH) {
					*ring_lagging_num = total_num;
					*ring_lagging_procs = total_procs;
				}
			}
			else {
				//��������
				send_procs_template(comm, &right_proc, RING_MULTICAST_LAGGING_NUM, RING_MULTICAST_LAGGING_PROCS, total_num, total_procs);
			}
			if (right_proc.comm_stage == FINISH) break;
		}
		break;
	}
	case 2:
	{
		while (1) {
			recv_procs_template(comm, &left_proc, RING_GATHER_LAGGING_NUM, RING_GATHER_LAGGING_PROCS, &recv_num, &recv_procs);
			if (left_proc.comm_stage == FINISH) break;
		}

		merge_procs(local_num, recv_num, local_procs, recv_procs, &total_num, &total_procs);
		*ring_lagging_num = total_num;
		*ring_lagging_procs = total_procs;

		//�ڶ��֣��ಥ
		left_proc.comm_stage = READY; right_proc.comm_stage = READY;
		while (1)
		{
			send_procs_template(comm, &right_proc, RING_MULTICAST_LAGGING_NUM, RING_MULTICAST_LAGGING_PROCS, total_num, total_procs);
			if (right_proc.comm_stage == FINISH) { break; }
		}
		break;
	}
	default: 
		break; 
	}
	return GATHER_SUCCESS;
}