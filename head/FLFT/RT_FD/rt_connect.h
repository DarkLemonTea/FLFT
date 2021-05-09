#pragma once

#ifndef HEADER_FILE
#define HEADER_FILE
#include "../fl_head.h"
#endif

#include "../Ring_FD/ring_connect.h"
#include "rt_count.h"

//=============================================================================================================================
//构造重连进程的链表
//=============================================================================================================================

//寻找最近有效的通信对象
int find_rt_vaild_neighbor(
	int my_rank,
	int comm_size,
	int ring_num,
	Recon_LL *head,
	int *neighbor_rank
) {
	if (head->next == NULL) { return 0; }

	int res = 0;
	int dis = ring_num;
	Recon_LL *p;
	p = head->next;
	while (p != NULL) {
		if (p->proc.comm_stage == CONNECT_FINISH &&
			(my_rank - p->proc.rank + ring_num) % ring_num < dis) {
			dis = (my_rank - p->proc.rank + ring_num) % ring_num;
			*neighbor_rank = p->proc.rank;
			res = 1;
		}
		p = p->next;
	}
	return res;
}

//=============================================================================================================================
//通信模板
//环、树方向的连接使用不同的tag
//=============================================================================================================================

//树被动通信模板
void tree_passive_comm_template(
	MPI_Comm comm,
	Comm_proc *proc
) {
	switch ((*proc).comm_stage)
	{
	case CONNECT_STATE1: {
		MPI_Iprobe((*proc).rank, TREE_CONNECT1, comm, &(*proc).flag, &(*proc).status);
		if ((*proc).flag == 1) {
			MPI_Recv(NULL, 0, MPI_BYTE, (*proc).rank, TREE_CONNECT1, comm, MPI_STATUSES_IGNORE);
			MPI_Isend(&(*proc).send, 1, MPI_INT, (*proc).rank, TREE_CONNECT2, comm, &(*proc).req);
			(*proc).comm_stage = CONNECT_STATE2;
		}
		break;
	}
	case CONNECT_STATE2: {
		MPI_Test(&(*proc).req, &(*proc).flag, MPI_STATUSES_IGNORE);
		if ((*proc).flag == 1) {
			(*proc).comm_stage = CONNECT_STATE3;
		}
		break;
	}
	case CONNECT_STATE3: {
		MPI_Iprobe((*proc).rank, TREE_CONNECT3, comm, &(*proc).flag, &(*proc).status);
		if ((*proc).flag == 1) {
			MPI_Recv(&(*proc).recv, 1, MPI_INT, (*proc).rank, TREE_CONNECT3, comm, MPI_STATUSES_IGNORE);
			if ((*proc).recv == (*proc).send) {
				MPI_Isend(NULL, 0, MPI_BYTE, (*proc).rank, TREE_CONNECT4, comm, &(*proc).req);
				(*proc).comm_stage = CONNECT_STATE4;
			}
			else if ((*proc).recv > (*proc).send) {
				(*proc).comm_stage = CONNECT_LAGGING;
			}
		}
		break;
	}
	case CONNECT_STATE4: {
		MPI_Test(&(*proc).req, &(*proc).flag, MPI_STATUSES_IGNORE);
		if ((*proc).flag == 1) {
			(*proc).comm_stage = CONNECT_FINISH;
		}
		break;
	}
	default:
		break;
	}
}

//树主动通信模板
void tree_active_comm_template(
	MPI_Comm comm,
	Comm_proc *proc
) {
	switch ((*proc).comm_stage)
	{
	case CONNECT_STATE1: {
		MPI_Isend(NULL, 0, MPI_BYTE, (*proc).rank, TREE_CONNECT1, comm, &(*proc).req);
		(*proc).comm_stage = CONNECT_STATE2;
		break;
	}
	case CONNECT_STATE2: {
		MPI_Test(&(*proc).req, &(*proc).flag, MPI_STATUSES_IGNORE);
		if ((*proc).flag == 1) {
			(*proc).comm_stage = CONNECT_STATE3;
		}
		break;
	}
	case CONNECT_STATE3: {
		MPI_Iprobe((*proc).rank, TREE_CONNECT2, comm, &(*proc).flag, &(*proc).status);
		if ((*proc).flag == 1) {
			MPI_Recv(&(*proc).recv, 1, MPI_INT, (*proc).rank, TREE_CONNECT2, comm, MPI_STATUSES_IGNORE);
			if ((*proc).recv == (*proc).send) {
				MPI_Isend(&(*proc).send, 1, MPI_INT, (*proc).rank, TREE_CONNECT3, comm, &(*proc).req);
				(*proc).comm_stage = CONNECT_STATE4;
			}
			else if ((*proc).recv > (*proc).send) {
				(*proc).comm_stage = CONNECT_LAGGING;
			}
		}
		break;
	}
	case CONNECT_STATE4: {
		MPI_Test(&(*proc).req, &(*proc).flag, MPI_STATUSES_IGNORE);
		if ((*proc).flag == 1) { (*proc).comm_stage = CONNECT_STATE5; }
		break;
	}
	case CONNECT_STATE5: {
		MPI_Iprobe((*proc).rank, TREE_CONNECT4, comm, &(*proc).flag, &(*proc).status);
		if ((*proc).flag == 1) {
			MPI_Recv(NULL, 0, MPI_BYTE, (*proc).rank, TREE_CONNECT4, comm, MPI_STATUSES_IGNORE);
			(*proc).comm_stage = CONNECT_FINISH;
		}
		break;
	}
	default:
		break;
	}
}

//=============================================================================================================================
//连接函数-new
//1—连接相邻有效对象
//2—计数
//3—超时重连
//4—解决错误连接
//=============================================================================================================================

//判断是否满足计数条件
void rt_ring_counting_demand(
	int my_rank,
	int comm_size,
	Ring_Tree rt,
	Recon_procs rp,
	Ring_Counter *c
) {
	//相邻进程都成功连接
	if (rt.ring.left_proc.comm_stage == CONNECT_FINISH &&
		rt.ring.right_proc.comm_stage == CONNECT_FINISH) {
		(*c).demand = 1;
		(*c).stage = STAGE1;
		(*c).r.left_proc.rank = rt.ring.left_proc.rank;
		(*c).r.right_proc.rank = rt.ring.right_proc.rank;
	}
	else if (rt.ring.left_proc.comm_stage != CONNECT_FINISH &&
		rt.ring.right_proc.comm_stage == CONNECT_FINISH) {
		(*c).demand = find_rt_vaild_neighbor(my_rank, comm_size, rt.ring_num,
			rp.left_head, &((*c).r.left_proc.rank));
		if ((*c).demand == 1) {
			(*c).r.right_proc.rank = rt.ring.right_proc.rank;
			(*c).stage = STAGE1;
		}
	}
	else if (rt.ring.left_proc.comm_stage == CONNECT_FINISH && 
		rt.ring.right_proc.comm_stage != CONNECT_FINISH) {
		(*c).demand = find_rt_vaild_neighbor(my_rank, comm_size, rt.ring_num,
			rp.right_head, &((*c).r.right_proc.rank));
		if ((*c).demand == 1) {
			(*c).r.left_proc.rank = rt.ring.left_proc.rank;
			(*c).stage = STAGE1;
		}
	}
	else {
		(*c).demand = find_rt_vaild_neighbor(my_rank, comm_size, rt.ring_num,
			rp.left_head, &((*c).r.left_proc.rank)) &&
			find_rt_vaild_neighbor(my_rank, comm_size, rt.ring_num, 
				rp.right_head, &((*c).r.right_proc.rank));
		if ((*c).demand == 1) {
			(*c).stage = STAGE1;
		}
	}
}

void rt_counting_demand(Tree t, Last_con_sta *lcs,RT_Counter *cou) {
	//对树的状态进行一次保存和截断。
	if (t.parent.comm_stage == CONNECT_FINISH) { 
		(*lcs).p_lcs = 1; 
	}
	if (t.left_child.comm_stage == CONNECT_FINISH) { 
		(*lcs).l_lcs = 1;
	}
	if (t.right_child.comm_stage == CONNECT_FINISH) { 
		(*lcs).r_lcs = 1;
	}

	(*cou).stage = STAGE1;
	(*cou).demand = 1;
}

int ring_tree_procs_connect(
	int detector_stage,
	FD_var fd_var,
	MPI_Comm comm,
	P_set ring_lagging_set,
	Ring_Tree *rt //用于存储当前有效连接对象
) {
	int my_rank, comm_size;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &comm_size);

	int min_arr_num = fd_var.release_rate * comm_size;
	if (min_arr_num > comm_size - 1) {
		min_arr_num = comm_size - 1;
	}
	else if (min_arr_num < comm_size - fd_var.max_fail_num) {
		min_arr_num = comm_size - fd_var.max_fail_num;
	}

	Recon_procs rp;
	init_recon_procs(&rp);
	RT_Counter cou;
	init_rt_counter((*rt).tree, &cou);
	Comm_proc ver_proc;

	int i, result = -1;
	int T_wait_out = 0;
	struct timeval start, end;
	double cost_time;
	gettimeofday(&start, NULL);
	rp.start.sec = start.tv_sec;
	rp.start.usec = start.tv_usec;
	
	rp.recently_tried_proc = (*rt).ring.right_proc.rank;
	
	init_lcs(&(*rt).lcs);

	

	while (1)
	{
		//环连接
		ring_passive_comm_template(comm, &(*rt).ring.left_proc);
		ring_active_comm_template(comm, &(*rt).ring.right_proc);

		//树连接
		tree_active_comm_template(comm, &(*rt).tree.parent);
		tree_passive_comm_template(comm, &(*rt).tree.left_child);
		tree_passive_comm_template(comm, &(*rt).tree.right_child);

		//判断是否重连
		if (T_wait_out == 0) {
			gettimeofday(&end, NULL);
			cost_time = 1000.0 * (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000.0;
			if (cost_time > fd_var.T_wait) {
				T_wait_out = 1;
				if ((*rt).ring.left_proc.comm_stage == CONNECT_FINISH) { rp.need_to_recon_l = 0; }
				else { rp.need_to_recon_l = 1; }

				if ((*rt).ring.right_proc.comm_stage == CONNECT_FINISH) { rp.need_to_recon_r = 0; }
				else { rp.need_to_recon_r = 1; }
			}
			
		}

		//重连
		if (T_wait_out == 1) {
			//左进程重连
			if (rp.need_to_recon_l == 1) {
				ring_passive_reconnect(detector_stage, comm, rp.left_head);
			}
			//右进程重连，如果超时，增加寻找下一个有效进程。
			if (rp.need_to_recon_r == 1) {
				gettimeofday(&end, NULL);
				cost_time = 1000.0 * (end.tv_sec - rp.start.sec) + (end.tv_usec - rp.start.usec) / 1000.0;
				if (cost_time > fd_var.T_wait) {
					//记录时间
					rp.start.sec = end.tv_sec;
					rp.start.usec = end.tv_usec;
					//找到下一个有效的通信对象
					rp.recently_tried_proc = find_rt_valid(comm_size, rp.recently_tried_proc, 
						(*rt).ring_num, ring_lagging_set, 'r');
					//将该进程加入链表中
					creat_new_recon_proc(rp.right_head, rp.recently_tried_proc, detector_stage);
				}
				ring_active_reconnect(comm, rp.right_head);
			}
			//如果直连进程到达，就不需要再重连，关闭重连通道
			if ((*rt).ring.left_proc.comm_stage == CONNECT_FINISH) { rp.need_to_recon_l = 0; }
			if ((*rt).ring.right_proc.comm_stage == CONNECT_FINISH) { rp.need_to_recon_r = 0; }
		}

		//判断参与计数的进程
		if (cou.demand == 0) { 		
			if (cou.rc.demand == 0) {
				rt_ring_counting_demand(my_rank, comm_size, (*rt), rp, &cou.rc);
				
				/*if (cou.rc.demand == 1) {
					printf("cou rank %d,ring left is %d,ring right is %d\n",
						my_rank, cou.rc.r.left_proc.rank, cou.rc.r.right_proc.rank);
				}*/
			}
			if (cou.rc.demand == 1) {
				if ((*rt).tree.parent.comm_stage == CONNECT_FINISH &&
					(*rt).tree.left_child.comm_stage == CONNECT_FINISH &&
					(*rt).tree.right_child.comm_stage == CONNECT_FINISH) {
					//如果树方向全部连接成功，初始化计数变量。
					//printf("rank %d tree connect success\n", my_rank);
					rt_counting_demand((*rt).tree, &(*rt).lcs, &cou);
				}
				else {
					gettimeofday(&end, NULL);
					cost_time = 1000.0 * (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000.0;
					if (cost_time > fd_var.T_wait) {
						//等待超过一定时间才允许计数
						rt_counting_demand((*rt).tree, &(*rt).lcs, &cou);
					}
				}
			}
		}
			
		//有计数需求，执行计数
		if (cou.demand == 1) { 
			ring_tree_arrived_procs_count(comm, my_rank, (*rt).lcs, &cou);
		}
		//如果计数完成
		if (cou.stage == FINISH) {
			/*if (my_rank == 0) {
				printf("rank %d count num is %d\n", my_rank, cou.sum);
			}*/
			
			if (cou.sum < min_arr_num) {
				//计数未达标，初始化counter，重新计数
				init_rt_counter((*rt).tree,&cou);
			}
			else {
				//达到计数条件，通过
				(*rt).ring.left_proc.rank = cou.rc.r.left_proc.rank;
				(*rt).ring.right_proc.rank = cou.rc.r.right_proc.rank;

				if (cou.sum == comm_size) {
					result = RT_CONNECT_SUCCESS;
				}
				else {

					result = RT_CONNECT_THROUGH;
				}
				break;
			}
		}

		//判断是否迟到
		if (is_lagging((*rt).ring, rp)) {
			//printf("rank %d res is connect lagging\n", my_rank);
			result = RT_CONNECT_LAGGING;
		}

		//判断是否有错误重连
		MPI_Iprobe(MPI_ANY_SOURCE, RING_MISCONNECTION, comm, &ver_proc.flag, &ver_proc.status);
		if (ver_proc.flag == 1) {
			//printf("rank %d misconnect\n", my_rank);
			MPI_Recv(NULL, 0, MPI_BYTE, ver_proc.status.MPI_SOURCE, 
				RING_MISCONNECTION, comm, MPI_STATUSES_IGNORE);
			result = RT_CONNECT_LAGGING;
		}
		//如果迟到，且右边进程有正确连接，发送通知
		if (result == RT_CONNECT_LAGGING) {
			//暂时取消错误重连的判断
			/*if (cou.rc.r.right_proc.rank > -1) {
				MPI_Isend(NULL, 0, MPI_BYTE, cou.rc.r.right_proc.rank,
					RING_MISCONNECTION, comm, &cou.rc.r.right_proc.req);
			}*/
			break;
		}

		gettimeofday(&end, NULL);
		cost_time = 1000.0 * (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000.0;
		//如果超过最大时限，返回失败
		if (cost_time > fd_var.T_max) {
			result = RT_CONNECT_FAILURE;
			break;
		}
	}

	//删除链表
	clear_recon_list(rp.left_head);
	clear_recon_list(rp.right_head);

	return result;
}
