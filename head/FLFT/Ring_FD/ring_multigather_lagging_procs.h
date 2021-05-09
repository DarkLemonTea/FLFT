#pragma once

#ifndef HEADER_FILE
#define HEADER_FILE
#include "ring_head.h"
#endif

//���ھۼ�lagging����
int multigather_lagging_procs(
	MPI_Comm comm,
	Ring ring,
	//output
	P_set *rlp
) {
	int my_rank, comm_size;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &comm_size);

	int ide, i;
	if (ring.left_proc.rank > my_rank) { ide = 0; }//��ʼ����
	else if (my_rank > ring.right_proc.rank) { ide = 2; }//ĩβ����
	else { ide = 1; }

	P_set recv_set, local_set, send_set;
	init_p_set(&recv_set);
	init_p_set(&send_set);
	init_p_set(&local_set);

	ring_statistics_lagging_procs(ring, my_rank, comm_size, &local_set);

	Comm_proc left_proc, right_proc, ver_proc;
	init_proc(ring.left_proc.rank, &left_proc);  
	init_proc(ring.right_proc.rank, &right_proc);

	switch (ide)
	{
	case 0:
	{	//��һ�ֻ�״�ಥ��������ʵ��	
		while (1) {
			send_procs_template(comm, &right_proc, 
				RING_GATHER_LAGGING_NUM, RING_GATHER_LAGGING_PROCS, local_set);
			if (right_proc.comm_stage == FINISH) break;
		}

		//�ڶ��ְ��ջ���
		left_proc.comm_stage = READY; right_proc.comm_stage = READY;
		while (1)
		{
			//��������
			if (left_proc.comm_stage != FINISH) {
				recv_procs_template(comm, &left_proc, 
					RING_MULTICAST_LAGGING_NUM, RING_MULTICAST_LAGGING_PROCS, rlp);
			}
			else {
				//��������
				send_procs_template(comm, &right_proc, 
					RING_MULTICAST_LAGGING_NUM, RING_MULTICAST_LAGGING_PROCS, *rlp);
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
				recv_procs_template(comm, &left_proc, 
					RING_GATHER_LAGGING_NUM, RING_GATHER_LAGGING_PROCS, &recv_set);
			}
			else {
				if (right_proc.comm_stage == READY) {
					merge_procs(local_set, recv_set, &send_set);
					send_procs_template(comm, &right_proc, 
						RING_GATHER_LAGGING_NUM, RING_GATHER_LAGGING_PROCS, send_set);
				}
				else {
					send_procs_template(comm, &right_proc, 
						RING_GATHER_LAGGING_NUM, RING_GATHER_LAGGING_PROCS, send_set);
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
				recv_procs_template(comm, &left_proc, 
					RING_MULTICAST_LAGGING_NUM, RING_MULTICAST_LAGGING_PROCS, rlp);
			}
			else {
				//��������
				send_procs_template(comm, &right_proc, 
					RING_MULTICAST_LAGGING_NUM, RING_MULTICAST_LAGGING_PROCS, *rlp);
			}
			if (right_proc.comm_stage == FINISH) break;
		}
		break;
	}
	case 2:
	{
		while (1) {
			recv_procs_template(comm, &left_proc, 
				RING_GATHER_LAGGING_NUM, RING_GATHER_LAGGING_PROCS, &recv_set);
			if (left_proc.comm_stage == FINISH) break;
		}

		merge_procs(local_set, recv_set, rlp);

		//�ڶ��֣��ಥ
		left_proc.comm_stage = READY; right_proc.comm_stage = READY;
		while (1)
		{
			send_procs_template(comm, &right_proc, 
				RING_MULTICAST_LAGGING_NUM, RING_MULTICAST_LAGGING_PROCS, *rlp);
			if (right_proc.comm_stage == FINISH) { break; }
		}
		break;
	}
	default: 
		break; 
	}

	return SUCCESS;
}