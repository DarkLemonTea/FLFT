#pragma once

#ifndef HEADER_FILE
#define HEADER_FILE
#include "rt_head.h"
#include "../Ring_FD/ring_count.h"
#endif

//=============================================================================================================================
//���������
//=============================================================================================================================
typedef struct rt_counter {
	int demand; //���ڼ�¼�Ƿ��м�������
	int stage;  //�ܵļ���ִ�н׶�

	int l_count; //���ӻ�����
	int r_count; //�Һ��ӻ�����
	int local_count; //����+�Һ���+����

	//�������ݽ���
	Tran_procs p_tp;	
	Tran_procs lc_tp;
	Tran_procs rc_tp;

	int sum; //�ܵ�����

	int ide; //�ж���������������Ҷ�ӻ�����֦
	Ring_Counter rc;
	Tree t;  //�������ݷ���
}RT_Counter;

void init_rt_counter(
	Tree t, 
	RT_Counter *c) {
	init_ring_counter(&(*c).rc);

	if (t.parent.rank == EMPTY) {
		(*c).ide = 0; //����
	}
	else if (t.left_child.rank == EMPTY && t.right_child.rank == EMPTY) {
		(*c).ide = 2; //��Ҷ
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
		//�׶�1���������׶�
		ring_arrived_procs_count(comm, my_rank, &(*cou).rc);
		
		if ((*cou).rc.stage == FINISH) {
			//��Ҫ��������Ĵ��ݱ������г�ʼ��
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
		//�׶�2����Ҷ�������ķ���reduce�ѵ����������
		switch ((*cou).ide)
		{
		case TREE_ROOT: {
			//�������������Һ��ӽ��̼���
			
			//���ӽ��̵�����,��ȡ����Ϊ����̣��������Ϊ�ҽ���
			count_transmit(comm, (*cou).t.left_child.rank, lcs.l_lcs, &(*cou).lc_tp,
				MISSING_REDUCE_COUNT_REQUEST, TREE_COUNT_REDUCE, &(*cou).l_count);
			
			//�Һ��ӽ��̵�����,��ȡ����Ϊ�ҽ��̣��������Ϊ�����
			count_transmit(comm, (*cou).t.right_child.rank, lcs.r_lcs, &(*cou).rc_tp,
				MISSING_REDUCE_COUNT_REQUEST, TREE_COUNT_REDUCE, &(*cou).r_count);

			if ((*cou).lc_tp.stage == FINISH && (*cou).rc_tp.stage == FINISH) {
				//����reduce���
				(*cou).sum = (*cou).l_count + (*cou).r_count + (*cou).rc.sum;
				
				(*cou).t.left_child.comm_stage = READY;
				(*cou).t.right_child.comm_stage = READY;
				(*cou).stage = STAGE3;
			}
			break;
		}
		case TREE_BRANCH: {
			//֦���������Һ��ӵļ������ϴ�

			if (!((*cou).lc_tp.stage == FINISH && (*cou).rc_tp.stage == FINISH)) {
				
				//���ӽ��̵�����,��ȡ����Ϊ����̣��������Ϊ�ҽ���
				count_transmit(comm, (*cou).t.left_child.rank, lcs.l_lcs, &(*cou).lc_tp,
					MISSING_REDUCE_COUNT_REQUEST, TREE_COUNT_REDUCE, &(*cou).l_count);

				//�Һ��ӽ��̵�����,��ȡ����Ϊ�ҽ��̣��������Ϊ�����
				count_transmit(comm, (*cou).t.right_child.rank, lcs.r_lcs, &(*cou).rc_tp,
					MISSING_REDUCE_COUNT_REQUEST, TREE_COUNT_REDUCE, &(*cou).r_count);
				
				if ((*cou).lc_tp.stage == FINISH && (*cou).rc_tp.stage == FINISH) {
					//����+�Һ���+����
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
			//Ҷ�ӣ�ֱ�����ϴ�
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
		//�׶�3���ɸ���Ҷ�ӵķ���bcast�ѵ����������
		switch ((*cou).ide)
		{
		case TREE_ROOT: {
			//�����̣������ݷ��͸����ӽ��̼���
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
			//һ����̣�����˫�׽��̷��������ݣ����͸����ӽ���
			if ((*cou).p_tp.stage != FINISH) {
				//����˫�׽������ݣ���ȡ����Ϊ����̣��������Ϊ�ҽ���
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
				//����˫�׽������ݣ���ȡ����Ϊ����̣��������Ϊ�ҽ���
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


