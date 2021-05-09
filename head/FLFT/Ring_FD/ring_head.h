#pragma once

#include"../base_head.h"

//————————————————————————————————————————————————————————
//返回结果
#define RING_CONNECT_SUCCESS 0
#define RING_CONNECT_THROUGH 1
#define RING_CONNECT_FAILURE 2
#define RING_CONNECT_LAGGING 3

//————————————————————————————————————————————————————————
//连接和重连
#define RING_CONNECT1 10002
#define RING_CONNECT2 10003
#define RING_CONNECT3 10004
#define RING_CONNECT4 10005
#define RING_RECONNECT1 10006
#define RING_RECONNECT2 10007
#define RING_RECONNECT3 10008
#define RING_RECONNECT4 10009

//计数
#define RING_COUNT_REDUCE 11000
#define RING_COUNT_BCAST 11001
//特殊案例，完成计数达到放行条件，此时有进程到达并成功建立连接，由左进程发送给右进程
#define RING_COUNT_SKIP 11002 

//————————————————————————————————————————————————————————
//迟到进程处理
#define RING_GATHER_LAGGING_NUM 20020
#define RING_GATHER_LAGGING_PROCS 20021
#define RING_MULTICAST_LAGGING_NUM 20022
#define RING_MULTICAST_LAGGING_PROCS 20023

#define RING_GATHER_REVIVE_NUM 20030
#define RING_GATHER_REVIVE_PROCS 20031
#define RING_MULTICAST_REVIVE_NUM 20032
#define RING_MULTICAST_REVIVE_PROCS 20033
#define RING_REVIVE_DATA_LAGGING_NUM 20040
#define RING_REVIVE_DATA_LAGGING_PROCS 20041

#define RING_LAGGING 21000
#define RING_REVIVE 21010
#define RING_LAGGING_PROBE 21011
#define RING_LAGGING_RESPOND 21012
#define RING_LAGGING_CONFIRM 21013
#define RING_RETRIEVE 21014
#define RING_CONFIRM 21015
#define RING_MISCONNECTION 21016


//=============================================================================================================================
//Ring
//=============================================================================================================================

typedef struct ring {
	Comm_proc left_proc;
	Comm_proc right_proc;
}Ring;

int left_proc(int comm_size, int my_rank) { return (comm_size + my_rank - 1) % comm_size; }
int right_proc(int comm_size, int my_rank) { return (my_rank + 1) % comm_size; }

//寻找最近有效对象
int find_ring_valid(int comm_size, int rank, P_set l, char l_or_r) {
	int tmp_proc;
	if (l_or_r == 'l') {
		tmp_proc = left_proc(comm_size, rank);
		while (in(l, tmp_proc)) { tmp_proc = left_proc(comm_size, tmp_proc); }		
	}
	else if (l_or_r == 'r') {
		tmp_proc = right_proc(comm_size, rank);
		while (in(l, tmp_proc)) { tmp_proc = right_proc(comm_size, tmp_proc); }	
	}
	return tmp_proc;
}

//初始化环
void init_ring(
	int comm_size,
	int my_rank,
	int stage,
	P_set l,
	Ring *ring
) {
	(*ring).left_proc.rank = find_ring_valid(comm_size, my_rank, l, 'l');
	(*ring).right_proc.rank = find_ring_valid(comm_size, my_rank, l, 'r');
	(*ring).left_proc.comm_stage = CONNECT_STATE1;
	(*ring).right_proc.comm_stage = CONNECT_STATE1;
	(*ring).left_proc.send = stage;
	(*ring).right_proc.send = stage;
}

void ring_statistics_lagging_procs(
	//input
	Ring ring,
	int my_rank,
	int comm_size,
	//output
	P_set *local
) {
	int i, num;
	int *procs;

	if (ring.left_proc.rank > my_rank) {
		num = my_rank + comm_size - 1 - ring.left_proc.rank;
	}
	else {
		num = my_rank - ring.left_proc.rank - 1;
	}

	if (num == 0) {
		(*local).num = 0;
		(*local).procs = NULL;
	}
	else {
		(*local).num = num;
		procs = (int*)calloc(num, sizeof(int));
		for (i = 0; i < num; i++) {
			procs[i] = (ring.left_proc.rank + 1 + i) % comm_size;
		}
		(*local).procs = procs;
	}
}

//===========================================================================================================================

////获得最后的进程号
//int get_the_last(
//	int comm_size,
//	P_set l
//) {
//	int last_proc = comm_size - 1;
//	while (in(l, last_proc)) {
//		last_proc -= 1;
//	}
//	return last_proc;
//}

//根据进程号和L，获取坐标，用于通信拓扑
int get_comm_index(
	int my_rank,
	P_set l,  //lagging procs
	int comm_size
) {
	int i, ind = my_rank;
	for (i = 0; i < l.num; i++) {
		if (l.procs[i] < my_rank) {
			ind -= 1;
		}
		else {
			break;
		}
	}
	return ind;
}

//根据L和index还原进程号
int get_real_rank(
	int index,
	P_set l, //lagging procs
	int comm_size
) {
	//已排序
	int i = 0;
	int real_num = index;
	while (i < l.num)
	{
		if (l.procs[i] <= real_num) {
			real_num = real_num + 1;
			i += 1;
		}
		else
		{
			break;
		}
	}
	return real_num % comm_size;
}

////去重合并数组
//void merge_revive_procs(
//	//input
//	int local_num,
//	int *local_procs,
//	int recv_num,
//	int *recv_procs,
//	//output
//	int *send_num,
//	int **send_procs
//) {
//	//合并
//	int num;
//	int *procs;
//	merge_procs(local_num, recv_num, local_procs, recv_procs, &num, &procs);
//	if (num == 0) {
//		*send_num = 0;
//		*send_procs = NULL;
//	}
//	else if (num == 1) {
//		*send_num = num;
//		*send_procs = procs;
//	}
//	else {
//		//排序
//		QuickSort(procs, 0, num - 1);
//		//去重
//		int i, new_num = 1;
//		int *new_procs;
//		new_procs = (int*)calloc(num, sizeof(int));
//		new_procs[0] = procs[0];
//		for (i = 1; i < num; i++) {
//			if (procs[i] == new_procs[new_num - 1]) {
//				continue;
//			}
//			else {
//				new_procs[new_num] = procs[i];
//				new_num += 1;
//			}
//		}
//		*send_num = new_num;
//		*send_procs = new_procs;
//		free(procs);
//	}
//}
//
