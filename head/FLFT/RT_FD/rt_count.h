#pragma once

#ifndef HEADER_FILE
#define HEADER_FILE
#include "rt_head.h"
#include "../Ring_FD/ring_count.h"
#endif

//=============================================================================================================================
//树方向计数
//=============================================================================================================================
typedef struct rt_counter {
	int demand; //用于记录是否有计数需求
	int stage;  //总的计数执行阶段

	int l_count; //左孩子环计数
	int r_count; //右孩子环计数
	int local_count; //左孩子+右孩子+自身

	//用于数据接收
	Tran_procs p_tp;	
	Tran_procs lc_tp;
	Tran_procs rc_tp;

	int sum; //总到达数

	int ide; //判断自身属于树根、叶子还是树枝
	Ring_Counter rc;
	Tree t;  //用于数据发送
}RT_Counter;

void init_rt_counter(
	Tree t, 
	RT_Counter *c) {
	init_ring_counter(&(*c).rc);

	if (t.parent.rank == EMPTY) {
		(*c).ide = 0; //树根
	}
	else if (t.left_child.rank == EMPTY && t.right_child.rank == EMPTY) {
		(*c).ide = 2; //树叶
	}
	else {
		(*c).ide = 1;
	}

	(*c).t.parent.rank = t.parent.rank;
	(*c).t.left_child.rank = t.left_child.rank;
	(*c).t.right_child.rank = t.right_child.rank;

	(*c).t.parent.comm_stage = READY;
	(*c).t.left_child.comm_stage = READY;
	(*c).t.right_child.comm_stage = READY;

	(*c).demand = 0;
	(*c).stage = READY;
	
	(*c).sum = 0;
	(*c).l_count = 0;
	(*c).r_count = 0;
	(*c).local_count = 0;
}

void ring_tree_arrived_procs_count(
	MPI_Comm comm,
	int my_rank,
	Last_con_sta lcs,
	RT_Counter *cou
) {
	switch ((*cou).stage)
	{
	case STAGE1: {
		//阶段1，环计数阶段
		ring_arrived_procs_count(comm, my_rank, &(*cou).rc);
		
		if ((*cou).rc.stage == FINISH) {
			//需要对树方向的传递变量进行初始化
			init_transmit_procs((*cou).t.parent.rank,
				(*cou).rc.r.left_proc.rank, (*cou).rc.r.right_proc.rank, &(*cou).p_tp);

			init_transmit_procs((*cou).t.left_child.rank,
				(*cou).rc.r.left_proc.rank, (*cou).rc.r.right_proc.rank, &(*cou).lc_tp);

			init_transmit_procs((*cou).t.right_child.rank,
				(*cou).rc.r.right_proc.rank, (*cou).rc.r.left_proc.rank, &(*cou).rc_tp);

			(*cou).stage = STAGE2;
		}
		break;
	}
	case STAGE2: {
		//阶段2，由叶子往根的方向，reduce已到达进程数据
		switch ((*cou).ide)
		{
		case TREE_ROOT: {
			//树根，接收左右孩子进程即可
			
			//左孩子进程的数据,索取对象为左进程，给予对象为右进程
			count_transmit(comm, (*cou).t.left_child.rank, lcs.l_lcs, &(*cou).lc_tp,
				MISSING_REDUCE_COUNT_REQUEST, TREE_COUNT_REDUCE, &(*cou).l_count);
			
			//右孩子进程的数据,索取对象为右进程，给予对象为左进程
			count_transmit(comm, (*cou).t.right_child.rank, lcs.r_lcs, &(*cou).rc_tp,
				MISSING_REDUCE_COUNT_REQUEST, TREE_COUNT_REDUCE, &(*cou).r_count);

			if ((*cou).lc_tp.stage == FINISH && (*cou).rc_tp.stage == FINISH) {
				//总数reduce完毕
				(*cou).sum = (*cou).l_count + (*cou).r_count + (*cou).rc.sum;
				
				(*cou).t.left_child.comm_stage = READY;
				(*cou).t.right_child.comm_stage = READY;
				(*cou).stage = STAGE3;
			}
			break;
		}
		case TREE_BRANCH: {
			//枝，接收左右孩子的计数，上传

			if (!((*cou).lc_tp.stage == FINISH && (*cou).rc_tp.stage == FINISH)) {
				
				//左孩子进程的数据,索取对象为左进程，给予对象为右进程
				count_transmit(comm, (*cou).t.left_child.rank, lcs.l_lcs, &(*cou).lc_tp,
					MISSING_REDUCE_COUNT_REQUEST, TREE_COUNT_REDUCE, &(*cou).l_count);

				//右孩子进程的数据,索取对象为右进程，给予对象为左进程
				count_transmit(comm, (*cou).t.right_child.rank, lcs.r_lcs, &(*cou).rc_tp,
					MISSING_REDUCE_COUNT_REQUEST, TREE_COUNT_REDUCE, &(*cou).r_count);
				
				if ((*cou).lc_tp.stage == FINISH && (*cou).rc_tp.stage == FINISH) {
					//左孩子+右孩子+本地
					(*cou).local_count = (*cou).l_count + (*cou).r_count + (*cou).rc.sum;
				}
			}
			else{
				if (lcs.p_lcs == 1) {
					send_num_template(comm, &(*cou).t.parent, TREE_COUNT_REDUCE, (*cou).local_count);
					//printf("rank %d local count is %d\n", my_rank, (*cou).local_count);
					
					if ((*cou).t.parent.comm_stage == FINISH) {
						(*cou).t.left_child.comm_stage = READY;
						(*cou).t.right_child.comm_stage = READY;
						(*cou).t.parent.comm_stage = READY;
						(*cou).stage = STAGE3;
					}
				}
				else {
					(*cou).t.left_child.comm_stage = READY;
					(*cou).t.right_child.comm_stage = READY;
					(*cou).t.parent.comm_stage = READY;
					(*cou).stage = STAGE3;
				}
			}
			break;
		}
		case TREE_LEAF: {
			//叶子，直接网上传
			if (lcs.p_lcs == 1) {
				send_num_template(comm, &(*cou).t.parent, TREE_COUNT_REDUCE, (*cou).rc.sum);
				if ((*cou).t.parent.comm_stage == FINISH) {
					(*cou).t.parent.comm_stage = READY;
					(*cou).stage = STAGE3;
				}
			}
			else {
				(*cou).t.parent.comm_stage = READY;
				(*cou).stage = STAGE3;
			}
			break;
		}
		default:
			break;
		}
		break;
	}
	case STAGE3: {
		//阶段3，由根往叶子的方向，bcast已到达进程数据
		switch ((*cou).ide)
		{
		case TREE_ROOT: {
			//根进程，将数据发送给孩子进程即可
			if (lcs.l_lcs == 0) {
				(*cou).t.left_child.comm_stage = FINISH;
			}
			else {
				send_num_template(comm, &(*cou).t.left_child, TREE_COUNT_BCAST, (*cou).sum);
			}

			if (lcs.r_lcs == 0) {
				(*cou).t.right_child.comm_stage = FINISH;
			}
			else {
				send_num_template(comm, &(*cou).t.right_child, TREE_COUNT_BCAST, (*cou).sum);
			}

			if ((*cou).t.left_child.comm_stage == FINISH && (*cou).t.right_child.comm_stage == FINISH) {
				(*cou).stage = FINISH;
			}
			break;
		}
		case TREE_BRANCH: {
			//一般进程，接收双亲进程发来的数据，发送给孩子进程
			if ((*cou).p_tp.stage != FINISH) {
				//接收双亲进程数据，索取对象为左进程，给予对象为右进程
				count_transmit(comm, (*cou).t.parent.rank, lcs.p_lcs, &(*cou).p_tp,
					MISSING_BCAST_COUNT_REQUEST, TREE_COUNT_BCAST, &(*cou).sum);
			}
			else {
				if (lcs.l_lcs == 0) {
					(*cou).t.left_child.comm_stage = FINISH;
				}
				else {
					send_num_template(comm, &(*cou).t.left_child, TREE_COUNT_BCAST, (*cou).sum);
				}

				if (lcs.r_lcs == 0) {
					(*cou).t.right_child.comm_stage = FINISH;
				}
				else {
					send_num_template(comm, &(*cou).t.right_child, TREE_COUNT_BCAST, (*cou).sum);
				}

				if ((*cou).t.left_child.comm_stage == FINISH &&
					(*cou).t.right_child.comm_stage == FINISH) {
					(*cou).stage = FINISH;
				}
			}
			break;
		}
		case TREE_LEAF: {
			if ((*cou).p_tp.stage != FINISH) {
				//接收双亲进程数据，索取对象为左进程，给予对象为右进程
				count_transmit(comm, (*cou).t.parent.rank, lcs.p_lcs, &(*cou).p_tp,
					MISSING_BCAST_COUNT_REQUEST, TREE_COUNT_BCAST, &(*cou).sum);
			}
			else {
				(*cou).stage = FINISH;
			}
		}
		default:
			break;
		}
		break;
	}
	default:
		break;
	}
}


