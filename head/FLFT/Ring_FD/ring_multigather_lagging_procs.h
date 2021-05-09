#pragma once

#ifndef HEADER_FILE
#define HEADER_FILE
#include "ring_head.h"
#endif

//环内聚集lagging进程
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
	if (ring.left_proc.rank > my_rank) { ide = 0; }//起始进程
	else if (my_rank > ring.right_proc.rank) { ide = 2; }//末尾进程
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
	{	//第一轮环状多播，非阻塞实现	
		while (1) {
			send_procs_template(comm, &right_proc, 
				RING_GATHER_LAGGING_NUM, RING_GATHER_LAGGING_PROCS, local_set);
			if (right_proc.comm_stage == FINISH) break;
		}

		//第二轮按照环形
		left_proc.comm_stage = READY; right_proc.comm_stage = READY;
		while (1)
		{
			//接收数组
			if (left_proc.comm_stage != FINISH) {
				recv_procs_template(comm, &left_proc, 
					RING_MULTICAST_LAGGING_NUM, RING_MULTICAST_LAGGING_PROCS, rlp);
			}
			else {
				//发送数组
				send_procs_template(comm, &right_proc, 
					RING_MULTICAST_LAGGING_NUM, RING_MULTICAST_LAGGING_PROCS, *rlp);
			}
			if (right_proc.comm_stage == FINISH) break;
		}
		break;
	}
	case 1:
	{
		//第一轮	
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

		//第二轮按照环形
		left_proc.comm_stage = READY; right_proc.comm_stage = READY;
		while (1)
		{
			//接收数组
			if (left_proc.comm_stage != FINISH) {
				recv_procs_template(comm, &left_proc, 
					RING_MULTICAST_LAGGING_NUM, RING_MULTICAST_LAGGING_PROCS, rlp);
			}
			else {
				//发送数组
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

		//第二轮，多播
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