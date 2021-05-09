#pragma once

#ifndef HEADER_FILE
#define HEADER_FILE
#include "../fl_head.h"
#endif

//���ھۼ�lagging����
int ring_multigather_lagging_procs(
	MPI_Comm comm,
	Ring ring,
	int ring_scale,
	int ring_num,
	//output
	P_set *rls
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

	rt_statistics_lagged_procs(ring, my_rank, comm_size,
		ring_scale, ring_num, &local_set);

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
			if (right_proc.comm_stage == FINISH) {
				left_proc.comm_stage = READY; right_proc.comm_stage = READY;
				break;
			}
		}

		//�ڶ��ְ��ջ���		
		while (1)
		{
			//��������
			if (left_proc.comm_stage != FINISH) {
				recv_procs_template(comm, &left_proc, 
					RING_MULTICAST_LAGGING_NUM, RING_MULTICAST_LAGGING_PROCS, rls);
			}
			else {
				//��������
				send_procs_template(comm, &right_proc, 
					RING_MULTICAST_LAGGING_NUM, RING_MULTICAST_LAGGING_PROCS, *rls);
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
					RING_MULTICAST_LAGGING_NUM, RING_MULTICAST_LAGGING_PROCS, rls);
			}
			else {
				//��������
				send_procs_template(comm, &right_proc, 
					RING_MULTICAST_LAGGING_NUM, RING_MULTICAST_LAGGING_PROCS, *rls);
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

		merge_procs(local_set, recv_set, rls);

		//�ڶ��֣��ಥ
		left_proc.comm_stage = READY; right_proc.comm_stage = READY;
		while (1)
		{
			send_procs_template(comm, &right_proc, 
				RING_MULTICAST_LAGGING_NUM, RING_MULTICAST_LAGGING_PROCS, *rls);
			if (right_proc.comm_stage == FINISH) { break; }
		}
		break;
	}
	default:
		break;
	}
	return SUCCESS;
}

int tree_multigather_lagging_procs(
	MPI_Comm comm,
	Ring_Tree rt,
	P_set rls,
	Last_con_sta lcs,
	//output
	P_set *ls
) {
	int my_rank, comm_size;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &comm_size);

	if (rt.tree.parent.rank == EMPTY &&
		rt.tree.left_child.rank == EMPTY &&
		rt.tree.right_child.rank == EMPTY) {
		return 0;
	}

	//��Ϊ������ķ��ͽ��̱���
	Comm_proc parent_proc, lchild_proc, rchild_proc;
	init_proc(rt.tree.parent.rank, &parent_proc);
	init_proc(rt.tree.left_child.rank, &lchild_proc);
	init_proc(rt.tree.right_child.rank, &rchild_proc);

	Tran_procs p_tp, lc_tp, rc_tp;
	//���ӽ��̵�����,��ȡ����Ϊ����̣��������Ϊ�ҽ���
	init_transmit_procs(rt.tree.left_child.rank,
		rt.ring.left_proc.rank, rt.ring.right_proc.rank, &lc_tp);
	//�Һ��ӽ��̵�����,��ȡ����Ϊ�ҽ��̣��������Ϊ�����
	init_transmit_procs(rt.tree.right_child.rank,
		rt.ring.right_proc.rank, rt.ring.left_proc.rank, &rc_tp);
	//����˫�׽������ݣ���ȡ����Ϊ����̣��������Ϊ�ҽ���	
	init_transmit_procs(rt.tree.parent.rank,
		rt.ring.left_proc.rank, rt.ring.right_proc.rank, &p_tp);

	P_set lrecv_set, rrecv_set, send_set;
	init_p_set(&lrecv_set);
	init_p_set(&rrecv_set);
	init_p_set(&send_set);

	/*
		�ۼ���Ҷ������������ۼ�
		�ಥ��������Ҷ�ӷ���ಥ
		ȱʧ�������ҽ�����Ҫ
	*/
	switch (rt.ide)
	{
	//================================================================================================================================================
	case TREE_ROOT: {
		//printf("rank %d, lchild is %d, cs is %d, rchild is %d, cs is %d\n", 
		//	my_rank, rt.tree.left_child.rank, rt.tree.left_child.comm_stage,
		//	rt.tree.right_child.rank, rt.tree.right_child.comm_stage);

		//���ۼ�����һ�����������Һ��ӵ�lagged������Ϣ
		while (1)
		{
			//���ӽ��̵�����,��ȡ����Ϊ����̣��������Ϊ�ҽ���
			data_transmit(comm, rt.tree.left_child, lcs.l_lcs, &lc_tp,
				MISSING_GATHER_DATA_REQUEST, TREE_GATHER_LAGGING_NUM, TREE_GATHER_LAGGING_PROCS,
				&lrecv_set);

			//�Һ��ӽ��̵�����,��ȡ����Ϊ�ҽ��̣��������Ϊ�����
			data_transmit(comm, rt.tree.right_child, lcs.r_lcs, &rc_tp,
				MISSING_GATHER_DATA_REQUEST, TREE_GATHER_LAGGING_NUM, TREE_GATHER_LAGGING_PROCS,
				&rrecv_set);

			if (lc_tp.stage == FINISH && rc_tp.stage == FINISH) {
				lchild_proc.comm_stage = READY; 
				rchild_proc.comm_stage = READY;
				break;
			}
		}

		//���ۼ����ڶ������ϲ�,�õ��ܵ�L����
		merge_tri_procs(rls, lrecv_set, rrecv_set, ls);
				
		//��������������������������������������������������������������������������������������������������������������������
		//���ಥ����һ��,���͸����Һ��ӽ���	
		while (1)
		{
			if (rt.tree.left_child.rank == EMPTY || lcs.l_lcs == 0) {
				lchild_proc.comm_stage = FINISH;
			}
			else {
				send_procs_template(comm, &lchild_proc, 
					TREE_MULTICAST_LAGGING_NUM, TREE_MULTICAST_LAGGING_PROCS, *ls);
			}

			if (rt.tree.right_child.rank == EMPTY || lcs.r_lcs == 0) {
				rchild_proc.comm_stage = FINISH;
			}
			else {
				send_procs_template(comm, &rchild_proc, 
					TREE_MULTICAST_LAGGING_NUM, TREE_MULTICAST_LAGGING_PROCS, *ls);
			}

			if (lchild_proc.comm_stage == FINISH && rchild_proc.comm_stage == FINISH) {
				//printf("rank %d has sent lagging procs to children\n", my_rank);
				break;
			}
		}
		//printf("rank %d finish tree multicast\n", my_rank);

		break;
	}
	//================================================================================================================================================
	case TREE_BRANCH: {

		//���ۼ�����һ�����������Һ��ӵ�lagged������Ϣ
		while (1)
		{
			//���ӽ��̵�����,��ȡ����Ϊ����̣��������Ϊ�ҽ���
			data_transmit(comm, rt.tree.left_child, lcs.l_lcs, &lc_tp,
				MISSING_GATHER_DATA_REQUEST, TREE_GATHER_LAGGING_NUM, TREE_GATHER_LAGGING_PROCS,
				&lrecv_set);

			//�Һ��ӽ��̵�����,��ȡ����Ϊ�ҽ��̣��������Ϊ�����
			data_transmit(comm, rt.tree.right_child, lcs.l_lcs, &rc_tp,
				MISSING_GATHER_DATA_REQUEST, TREE_GATHER_LAGGING_NUM, TREE_GATHER_LAGGING_PROCS,
				&rrecv_set);

			if (lc_tp.stage == FINISH && rc_tp.stage == FINISH) {
				lchild_proc.comm_stage = READY;
				rchild_proc.comm_stage = READY;
				break;
			}
		}
		//printf("rank %d finish tree gather\n", my_rank);

		//���ۼ����ڶ������ϲ������ϴ���	
		merge_tri_procs(rls, lrecv_set, rrecv_set, &send_set);
		//printf("rank %d merge lagged num is %d\n", my_rank, (*ls).num);

		while (1)
		{
			if (rt.tree.parent.rank == EMPTY || lcs.p_lcs == 0) {
				break;
			}
			else {
				send_procs_template(comm, &parent_proc, 
					TREE_GATHER_LAGGING_NUM, TREE_GATHER_LAGGING_PROCS, send_set);
				if (parent_proc.comm_stage == FINISH) {
					parent_proc.comm_stage = READY;
					break;
				}
			}
		}
		//��������������������������������������������������������������������������������������������������������������������
		//���ಥ����һ��,��������˫�׽��̵�L����
		while (1)
		{
			//����˫�׽������ݣ���ȡ����Ϊ����̣��������Ϊ�ҽ���
			data_transmit(comm, rt.tree.parent, lcs.p_lcs, &p_tp,
				MISSING_MULTICAST_DATA_REQUEST, TREE_MULTICAST_LAGGING_NUM, TREE_MULTICAST_LAGGING_PROCS,
				ls);
			if (p_tp.stage == FINISH) { 
				lchild_proc.comm_stage = READY; rchild_proc.comm_stage = READY;
				break; 
			}
		}
		//printf("rank %d recv parent multicast\n", my_rank);

		//���ಥ���ڶ��������͸����Һ���
		
		while (1)
		{
			if (rt.tree.left_child.rank == EMPTY || lcs.l_lcs == 0) {
				lchild_proc.comm_stage = FINISH;
			}
			else {
				send_procs_template(comm, &lchild_proc, 
					TREE_MULTICAST_LAGGING_NUM, TREE_MULTICAST_LAGGING_PROCS, *ls);
			}

			if (rt.tree.right_child.rank == EMPTY || lcs.r_lcs == 0) {
				rchild_proc.comm_stage = FINISH;
			}
			else {
				send_procs_template(comm, &rchild_proc, 
					TREE_MULTICAST_LAGGING_NUM, TREE_MULTICAST_LAGGING_PROCS, *ls);
			}

			if (lchild_proc.comm_stage == FINISH && rchild_proc.comm_stage == FINISH) {
				break;
			}
		}
		//printf("rank %d finish tree multicast\n", my_rank);
		break;
	}
	//================================================================================================================================================
	case TREE_LEAF: {
		//���ۼ�����һ�������ϴ��ݻ���lagged��������
		while (1)
		{
			if (rt.tree.parent.rank == EMPTY || lcs.p_lcs == 0) {
				break;
			}
			else {
				send_procs_template(comm, &parent_proc, 
					TREE_GATHER_LAGGING_NUM, TREE_GATHER_LAGGING_PROCS, rls);
				if (parent_proc.comm_stage == FINISH) {
					//printf("rank %d has sent local lagging procs to parent\n", my_rank);
					parent_proc.comm_stage = READY;
					break;
				}
			}
		}
		//printf("rank %d finish tree gather\n", my_rank);

		//��������������������������������������������������������������������������������������������������������������������
		//���ಥ����һ��,��������˫�׽��̵�L���ϣ�δ���ӳɹ����ҽ�����ȡ
		
		//printf("rank %d, parent is %d, connect status is %d\n", my_rank, rt.tree.parent.rank, rt.tree.parent.comm_stage);
		while (1)
		{
			//����˫�׽������ݣ���ȡ����Ϊ����̣��������Ϊ�ҽ���		
			data_transmit(comm, rt.tree.parent, lcs.p_lcs, &p_tp, MISSING_MULTICAST_DATA_REQUEST, 
				TREE_MULTICAST_LAGGING_NUM, TREE_MULTICAST_LAGGING_PROCS, ls);
			if (p_tp.stage == FINISH) { break; }
		}
		//printf("rank %d finish tree multicast\n", my_rank);
		break;
	}
			//================================================================================================================================================	
	default:
		break;
	}

	//printf("rank %d tree lagging num is %d\n", my_rank, (*ls).num);
	return SUCCESS;
}

int rt_multigather_lagging_procs(
	MPI_Comm comm,
	Ring_Tree rt,
	//output
	Detector *d
) {	
	int my_rank;
	MPI_Comm_rank(comm, &my_rank);

	ring_multigather_lagging_procs(comm, rt.ring, rt.ring_scale, rt.ring_num, &(*d).Ring_lagging_procs);
	//printf("rank %d ring lagging num is %d\n", my_rank, (*d).Ring_lagging_procs.num);

	tree_multigather_lagging_procs(comm, rt, (*d).Ring_lagging_procs, (*d).rt.lcs, &(*d).Lagging_procs);
	//printf("rank %d tree lagging num is %d\n", my_rank, (*d).Lagging_procs.num);

	return SUCCESS;
}