#pragma once

#ifndef HEADER_FILE
#define HEADER_FILE
#include "../fl_head.h"
#endif

//环内聚集lagging进程
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
	if (ring.left_proc.rank > my_rank) { ide = 0; }//起始进程
	else if (my_rank > ring.right_proc.rank) { ide = 2; }//末尾进程
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
	{	//第一轮环状多播，非阻塞实现	
		while (1) {
			send_procs_template(comm, &right_proc, 
				RING_GATHER_LAGGING_NUM, RING_GATHER_LAGGING_PROCS, local_set);
			if (right_proc.comm_stage == FINISH) {
				left_proc.comm_stage = READY; right_proc.comm_stage = READY;
				break;
			}
		}

		//第二轮按照环形		
		while (1)
		{
			//接收数组
			if (left_proc.comm_stage != FINISH) {
				recv_procs_template(comm, &left_proc, 
					RING_MULTICAST_LAGGING_NUM, RING_MULTICAST_LAGGING_PROCS, rls);
			}
			else {
				//发送数组
				send_procs_template(comm, &right_proc, 
					RING_MULTICAST_LAGGING_NUM, RING_MULTICAST_LAGGING_PROCS, *rls);
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
					RING_MULTICAST_LAGGING_NUM, RING_MULTICAST_LAGGING_PROCS, rls);
			}
			else {
				//发送数组
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

		//第二轮，多播
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

	//作为树方向的发送进程变量
	Comm_proc parent_proc, lchild_proc, rchild_proc;
	init_proc(rt.tree.parent.rank, &parent_proc);
	init_proc(rt.tree.left_child.rank, &lchild_proc);
	init_proc(rt.tree.right_child.rank, &rchild_proc);

	Tran_procs p_tp, lc_tp, rc_tp;
	//左孩子进程的数据,索取对象为左进程，给予对象为右进程
	init_transmit_procs(rt.tree.left_child.rank,
		rt.ring.left_proc.rank, rt.ring.right_proc.rank, &lc_tp);
	//右孩子进程的数据,索取对象为右进程，给予对象为左进程
	init_transmit_procs(rt.tree.right_child.rank,
		rt.ring.right_proc.rank, rt.ring.left_proc.rank, &rc_tp);
	//接收双亲进程数据，索取对象为左进程，给予对象为右进程	
	init_transmit_procs(rt.tree.parent.rank,
		rt.ring.left_proc.rank, rt.ring.right_proc.rank, &p_tp);

	P_set lrecv_set, rrecv_set, send_set;
	init_p_set(&lrecv_set);
	init_p_set(&rrecv_set);
	init_p_set(&send_set);

	/*
		聚集由叶子往树根方向聚集
		多播由树根往叶子方向多播
		缺失数据向右进程索要
	*/
	switch (rt.ide)
	{
	//================================================================================================================================================
	case TREE_ROOT: {
		//printf("rank %d, lchild is %d, cs is %d, rchild is %d, cs is %d\n", 
		//	my_rank, rt.tree.left_child.rank, rt.tree.left_child.comm_stage,
		//	rt.tree.right_child.rank, rt.tree.right_child.comm_stage);

		//（聚集）第一步，接收左右孩子的lagged进程信息
		while (1)
		{
			//左孩子进程的数据,索取对象为左进程，给予对象为右进程
			data_transmit(comm, rt.tree.left_child, lcs.l_lcs, &lc_tp,
				MISSING_GATHER_DATA_REQUEST, TREE_GATHER_LAGGING_NUM, TREE_GATHER_LAGGING_PROCS,
				&lrecv_set);

			//右孩子进程的数据,索取对象为右进程，给予对象为左进程
			data_transmit(comm, rt.tree.right_child, lcs.r_lcs, &rc_tp,
				MISSING_GATHER_DATA_REQUEST, TREE_GATHER_LAGGING_NUM, TREE_GATHER_LAGGING_PROCS,
				&rrecv_set);

			if (lc_tp.stage == FINISH && rc_tp.stage == FINISH) {
				lchild_proc.comm_stage = READY; 
				rchild_proc.comm_stage = READY;
				break;
			}
		}

		//（聚集）第二步，合并,得到总的L集合
		merge_tri_procs(rls, lrecv_set, rrecv_set, ls);
				
		//——————————————————————————————————————————————————————————
		//（多播）第一步,发送给左右孩子进程	
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

		//（聚集）第一步，接收左右孩子的lagged进程信息
		while (1)
		{
			//左孩子进程的数据,索取对象为左进程，给予对象为右进程
			data_transmit(comm, rt.tree.left_child, lcs.l_lcs, &lc_tp,
				MISSING_GATHER_DATA_REQUEST, TREE_GATHER_LAGGING_NUM, TREE_GATHER_LAGGING_PROCS,
				&lrecv_set);

			//右孩子进程的数据,索取对象为右进程，给予对象为左进程
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

		//（聚集）第二步，合并后，向上传递	
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
		//——————————————————————————————————————————————————————————
		//（多播）第一步,接收来自双亲进程的L集合
		while (1)
		{
			//接收双亲进程数据，索取对象为左进程，给予对象为右进程
			data_transmit(comm, rt.tree.parent, lcs.p_lcs, &p_tp,
				MISSING_MULTICAST_DATA_REQUEST, TREE_MULTICAST_LAGGING_NUM, TREE_MULTICAST_LAGGING_PROCS,
				ls);
			if (p_tp.stage == FINISH) { 
				lchild_proc.comm_stage = READY; rchild_proc.comm_stage = READY;
				break; 
			}
		}
		//printf("rank %d recv parent multicast\n", my_rank);

		//（多播）第二步，发送给左右孩子
		
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
		//（聚集）第一步，向上传递环内lagged进程数据
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

		//——————————————————————————————————————————————————————————
		//（多播）第一步,接收来自双亲进程的L集合，未连接成功，右进程索取
		
		//printf("rank %d, parent is %d, connect status is %d\n", my_rank, rt.tree.parent.rank, rt.tree.parent.comm_stage);
		while (1)
		{
			//接收双亲进程数据，索取对象为左进程，给予对象为右进程		
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