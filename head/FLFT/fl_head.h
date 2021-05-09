#pragma once

#ifndef HEADER_FILE
#define HEADER_FILE
#include "../FLFT/Ring_FD/ring_head.h"
#endif

#include "../FLFT/RT_FD/rt_head.h"

//----------------------------------------------------------------------------------------------
#define FD_SUCCESS 0
#define FD_FAILURE 1
#define FD_REVIVE 2
#define FD_PASS 3

#define GATHER_SUCCESS 0
#define CONFIRM_LAGGING 1
#define UNKNOWN_ERROR 2

#define TYPE_NONE 9			//不容错
#define TYPE_RING 10		//环拓扑
#define TYPE_RING_TREE 11		//环-树拓扑
#define TYPE_RING_BUTTERFLY 12		//环-蝶形拓扑

//等待时间与放行比例
typedef struct FD_variable {
	double release_rate;
	int max_fail_num;
	double T_max;
	double T_wait;
	double T_retrieve;
}FD_var;

void init_fd_var(
	double release_rate, //放行率
	int max_fail_num, //最大故障进程数量
	double T_max,        //wall time 超过该时间，返回错误
	double T_wait,       //wait time 超过该时间，并且达到放行率，完成故障检测
	double T_retrieve,   //进程捡回的等待时间
	FD_var *fd
) {
	//输入单位是秒，转化为毫秒
	(*fd).T_max = T_max * 1000;
	(*fd).T_wait = T_wait * 1000;
	(*fd).T_retrieve = T_retrieve * 1000;
	(*fd).release_rate = release_rate;
	(*fd).max_fail_num = max_fail_num;
}

//=============================================================================================================================
//通用的detector，根据通信结构不同选择初始化方式
//=============================================================================================================================

typedef struct Set_pointers {
	/*通用结构*/
	int type;					 //拓扑类型
	int current_stage;			 //当前阶段
	P_set Lagging_procs;         //迟到进程列表
	Revive_LL *local_re_phead;   //存储恢复进程的临时链表	
	P_set Revive_procs;          //恢复进程列表

	P_set Ring_lagging_procs;    //同环内迟到进程列表，用于快速初始化
	P_set Ring_Revive_procs;     //同环内恢复进程列表

	//环
	Ring ring;
	
	//环-树
	Ring_Tree rt;
	
}Detector;

void init_detector(
	int type,
	MPI_Comm comm,
	int d_stage,
	int ring_scale,  //if type is ring, ring_scale can be anything
	Detector *sp
) {
	int comm_size;
	int my_rank;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &comm_size);

	//初始化
	(*sp).current_stage = d_stage;
	(*sp).type = type;

	init_p_set(&(*sp).Lagging_procs);
	init_p_set(&(*sp).Revive_procs);

	//创建一个头结点，不代表任何数据
	(*sp).local_re_phead = (Revive_LL*)malloc(sizeof(Revive_LL));
	(*sp).local_re_phead->revive_proc.rank = -1;
	(*sp).local_re_phead->next = NULL;

	switch (type)
	{
	case TYPE_RING: {
		init_ring(comm_size, my_rank, d_stage, (*sp).Lagging_procs, &(*sp).ring);
		break;
	}
	case TYPE_RING_TREE: {
		init_p_set(&(*sp).Ring_lagging_procs);
		init_p_set(&(*sp).Ring_Revive_procs);

		init_ring_tree(comm_size, my_rank, ring_scale, d_stage,
			(*sp).Ring_lagging_procs, &(*sp).rt);

		if ((*sp).rt.tree.parent.rank == EMPTY) {
			(*sp).rt.ide = 0;
		}
		else if ((*sp).rt.tree.left_child.rank == EMPTY && (*sp).rt.tree.right_child.rank == EMPTY) {
			(*sp).rt.ide = 2;
		}
		else {
			(*sp).rt.ide = 1;
		}
		
		break;
	}
	default:
		break;
	}
}

//=============================================================================================================================
//存放revive进程的链表
//=============================================================================================================================

//带头结点（头结点不存储数据），头插法
void creat_and_init(Revive_LL *head, int rank, int stage) {
	Revive_LL *node;
	node = (Revive_LL*)malloc(sizeof(Revive_LL));
	node->revive_proc.rank = rank;
	node->revive_proc.send = stage;
	node->revive_proc.comm_stage = STAGE1;

	node->next = head->next;
	head->next = node;
}

//查找
int linked_list_in(Revive_LL *head, int target) {
	if (head->next == NULL) {
		return 0;
	}
	else {
		Revive_LL *p;
		p = head->next;
		while (p != NULL) {
			if (p->revive_proc.rank == target) {
				printf("rank %d in lagging procs\n", p->revive_proc.rank);
				return 1;
			}
			p = p->next;
		}
		return 0;
	}
}

//=============================================================================================================================
//清理消息
//=============================================================================================================================
void clear_message0(
	MPI_Comm comm,
	int tag,
	Ring ring
) {
	Comm_proc clear;
	while (1)
	{
		MPI_Iprobe(MPI_ANY_SOURCE, tag, comm, &clear.flag, &clear.status);
		if (clear.flag == 1 && clear.status.MPI_SOURCE != ring.left_proc.rank) {
			clear.rank = clear.status.MPI_SOURCE;
			MPI_Recv(NULL, 0, MPI_BYTE, clear.rank, tag, comm, MPI_STATUSES_IGNORE);
		}
		else {
			break;
		}
	}
}

void clear_message1(
	MPI_Comm comm,
	int tag
) {
	Comm_proc clear;
	while (1)
	{
		MPI_Iprobe(MPI_ANY_SOURCE, tag, comm, &clear.flag, &clear.status);
		if (clear.flag == 1) {
			clear.rank = clear.status.MPI_SOURCE;
			MPI_Recv(NULL, 0, MPI_BYTE, clear.rank, tag, comm, MPI_STATUSES_IGNORE);
		}
		else {
			break;
		}
	}
}

void clear_message2(
	MPI_Comm comm,
	int tag
) {
	Comm_proc clear;
	while (1)
	{
		MPI_Iprobe(MPI_ANY_SOURCE, tag, comm, &clear.flag, &clear.status);
		if (clear.flag == 1) {
			clear.rank = clear.status.MPI_SOURCE;
			MPI_Recv(&clear.recv, 1, MPI_INT, clear.rank, tag, comm, MPI_STATUSES_IGNORE);
		}
		else {
			break;
		}
	}
}

//===========================================================================================================================
//环的子操作
//===========================================================================================================================
//随机选一个同环进程
int get_subring_rand_proc(
	int my_rank,
	int ring_scale,
	int comm_size
) {
	int proc = rand() % comm_size;
	proc = proc % ring_scale;
	while (proc == my_rank ||
		proc == left_rt_proc(comm_size, my_rank, ring_scale) ||
		proc == right_rt_proc(comm_size, my_rank, ring_scale)) {
		proc = rand() % comm_size;
	}
	return proc;
}

//统计迟到进程（进程自身负责的部分）
void sub_ring_statistics_lagged_procs(
	//input
	Ring ring,
	int my_rank,
	int comm_size,
	int ring_scale,
	int ring_num,
	//output
	P_set *local
) {
	if (ring.left_proc.rank > my_rank) {
		(*local).num = my_rank + ring_num - 1 - ring.left_proc.rank;
	}
	else {
		(*local).num = my_rank - ring.left_proc.rank - 1;
	}

	if ((*local).num == 0) {
		(*local).num = 0;
		(*local).procs = NULL;
	}
	else {
		int i;
		int *procs;
		procs = (int*)calloc((*local).num, sizeof(int));
		for (i = 0; i < (*local).num; i++) {
			procs[i] = (my_rank / ring_scale) * ring_scale +
				(ring.left_proc.rank % ring_scale + 1 + i + ring_num) % ring_num;
		}
		(*local).procs = procs;
	}
}

//从lagging进程中获取同环内进程
void get_sub_ring_procs(
	int my_rank,
	int ring_scale,
	P_set p,
	P_set *rp
) {
	if (p.num == 0) {
		(*rp).num = 0;
		(*rp).procs = NULL;
		return;
	}

	(*rp).num = 0;
	int i = 0, ind = 0;
	for (i = 0; i < p.num; i++) {
		//同商同环
		if (my_rank / ring_scale == p.procs[i] / ring_scale) {
			(*rp).num += 1;
		}
	}
	(*rp).procs = (int*)calloc((*rp).num, sizeof(int));
	for (i = 0; i < p.num; i++) {
		//同商同环
		if (my_rank / ring_scale == p.procs[i] / ring_scale) {
			(*rp).procs[ind] = p.procs[i];
			ind += 1;
		}
	}
}

//===========================================================================================================================
//树、蝶方向数据传输、增补
//===========================================================================================================================
typedef struct transmit_procs {
	int stage;            //传输阶段
	Comm_proc dir_proc;  //树、环方向的对象，初始化的通信变量
	Comm_proc req_proc;   //索取对象，缺失时找该对象索要
	Comm_proc give_proc;  //给予对象，进程取得数据后，可以接受缺失请求
}Tran_procs;

void init_transmit_procs(
	int dir_rank,
	int req_rank,
	int give_rank,
	Tran_procs *tp
) {
	(*tp).stage = READY;

	(*tp).dir_proc.rank = dir_rank;
	(*tp).req_proc.rank = req_rank;
	(*tp).give_proc.rank = give_rank;

	(*tp).dir_proc.comm_stage = READY;
	(*tp).req_proc.comm_stage = READY;
	(*tp).give_proc.comm_stage = READY;
}

//数据传输
void data_transmit(
	MPI_Comm comm,
	Comm_proc dir_proc,
	int con_status,
	Tran_procs *tp,
	int missing_tag,
	int num_tag,
	int procs_tag,
	P_set *r
) {
	switch ((*tp).stage) {
	case READY: {
		(*tp).dir_proc.comm_stage = READY;
		(*tp).req_proc.comm_stage = READY;
		(*tp).give_proc.comm_stage = READY;

		if (dir_proc.rank == EMPTY) {
			(*tp).stage = FINISH;
		}
		else if (con_status == 1) {
			(*tp).req_proc.send = 0; //不需要填补缺失数据
			(*tp).stage = RECV_STAGE;
		}
		else if (con_status == 0) {
			(*tp).req_proc.send = 1; //表示有填补缺失数据请求
			(*tp).stage = REQ_CONFIRM_STAGE;
		}
		break;
	}
	case RECV_STAGE: {
		recv_procs_template(comm, &(*tp).dir_proc, num_tag, procs_tag, r);
		if ((*tp).dir_proc.comm_stage == FINISH) {
			(*tp).stage = REQ_CONFIRM_STAGE;
		}
		break;
	}
	case REQ_CONFIRM_STAGE: {
		//向索取放发送信号
		switch ((*tp).req_proc.comm_stage) {
		case READY: {
			MPI_Isend(&(*tp).req_proc.send, 1, MPI_INT, (*tp).req_proc.rank,
				missing_tag, comm, &(*tp).req_proc.req);
			(*tp).req_proc.comm_stage = STAGE1;
			break;
		}
		case STAGE1: {
			MPI_Test(&(*tp).req_proc.req, &(*tp).req_proc.flag, MPI_STATUSES_IGNORE);
			if ((*tp).req_proc.flag == 1) {
				(*tp).req_proc.comm_stage = FINISH;
			}
			break;
		}
		default:
			break;
		}

		//向接收给予方信号
		switch ((*tp).give_proc.comm_stage) {
		case READY: {
			MPI_Iprobe((*tp).give_proc.rank, missing_tag, comm, &(*tp).give_proc.flag, MPI_STATUSES_IGNORE);
			if ((*tp).give_proc.flag == 1) {
				MPI_Recv(&(*tp).give_proc.recv, 1, MPI_INT, (*tp).give_proc.rank,
					missing_tag, comm, MPI_STATUSES_IGNORE);
				(*tp).give_proc.comm_stage = FINISH;
			}
			break;
		}
		default:
			break;
		}

		if ((*tp).req_proc.comm_stage == FINISH && (*tp).give_proc.comm_stage == FINISH) {
			(*tp).req_proc.comm_stage = READY;
			(*tp).give_proc.comm_stage = READY;
			if ((*tp).req_proc.send == 1) {
				(*tp).stage = REQ_RECV_STAGE;
			}
			else if ((*tp).give_proc.recv == 1) {
				(*tp).stage = GIVE_STAGE;
			}
			else {
				(*tp).stage = FINISH;
			}
		}
		break;
	}
	case REQ_RECV_STAGE: {
		recv_procs_template(comm, &(*tp).req_proc, num_tag, procs_tag, r);
		if ((*tp).req_proc.comm_stage == FINISH) {
			if ((*tp).give_proc.recv == 1) {
				(*tp).give_proc.comm_stage = READY;
				(*tp).stage = GIVE_STAGE;
			}
			else {
				(*tp).stage = FINISH;
			}
		}
		break;
	}
	case GIVE_STAGE: {
		send_procs_template(comm, &(*tp).give_proc, num_tag, procs_tag, *r);
		if ((*tp).give_proc.comm_stage == FINISH) {
			(*tp).stage = FINISH;
		}
		break;
	}
	default:
		break;
	}
}

//计数传输
void count_transmit(
	MPI_Comm comm,
	int dir_rank,
	int con_status,
	Tran_procs *tp,
	int missing_tag,
	int num_tag,
	int *recv_num
) {
	switch ((*tp).stage) {
	case READY: {
		(*tp).dir_proc.comm_stage = READY;
		(*tp).req_proc.comm_stage = READY;
		(*tp).give_proc.comm_stage = READY;

		if (dir_rank == EMPTY) {
			(*tp).stage = FINISH;
		}
		else if (con_status == 1) {
			(*tp).req_proc.send = 0; //不需要填补缺失数据
			(*tp).stage = RECV_STAGE;
		}
		else if (con_status == 0) {
			(*tp).req_proc.send = 1; //表示有填补缺失数据请求
			(*tp).stage = REQ_CONFIRM_STAGE;
		}
		break;
	}
	case RECV_STAGE: {
		MPI_Iprobe((*tp).dir_proc.rank, num_tag, comm, &(*tp).dir_proc.flag, MPI_STATUSES_IGNORE);
		if ((*tp).dir_proc.flag == 1) {
			MPI_Recv(recv_num, 1, MPI_INT, (*tp).dir_proc.rank, num_tag, comm, MPI_STATUSES_IGNORE);
			(*tp).stage = REQ_CONFIRM_STAGE;
		}
		break;
	}
	case REQ_CONFIRM_STAGE: {
		//向索取放发送信号
		switch ((*tp).req_proc.comm_stage) {
		case READY: {
			MPI_Isend(&(*tp).req_proc.send, 1, MPI_INT, (*tp).req_proc.rank,
				missing_tag, comm, &(*tp).req_proc.req);
			(*tp).req_proc.comm_stage = STAGE1;
			break;
		}
		case STAGE1: {
			MPI_Test(&(*tp).req_proc.req, &(*tp).req_proc.flag, MPI_STATUSES_IGNORE);
			if ((*tp).req_proc.flag == 1) {
				(*tp).req_proc.comm_stage = FINISH;
			}
			break;
		}
		default:
			break;
		}

		//向接收给予方信号
		switch ((*tp).give_proc.comm_stage) {
		case READY: {
			MPI_Iprobe((*tp).give_proc.rank, missing_tag, comm, &(*tp).give_proc.flag, MPI_STATUSES_IGNORE);
			if ((*tp).give_proc.flag == 1) {
				MPI_Recv(&(*tp).give_proc.recv, 1, MPI_INT, (*tp).give_proc.rank,
					missing_tag, comm, MPI_STATUSES_IGNORE);
				(*tp).give_proc.comm_stage = FINISH;
			}
			break;
		}
		default:
			break;
		}

		if ((*tp).req_proc.comm_stage == FINISH && (*tp).give_proc.comm_stage == FINISH) {
			(*tp).req_proc.comm_stage = READY;
			(*tp).give_proc.comm_stage = READY;
			if ((*tp).req_proc.send == 1) {
				(*tp).stage = REQ_RECV_STAGE;
			}
			else if ((*tp).give_proc.recv == 1) {
				(*tp).stage = GIVE_STAGE;
			}
			else {
				(*tp).stage = FINISH;
			}
		}
		break;
	}
	case REQ_RECV_STAGE: {
		MPI_Iprobe((*tp).req_proc.rank, num_tag, comm, &(*tp).req_proc.flag, MPI_STATUSES_IGNORE);
		if ((*tp).req_proc.flag == 1) {
			MPI_Recv(recv_num, 1, MPI_INT, (*tp).req_proc.rank, num_tag, comm, MPI_STATUSES_IGNORE);
			if ((*tp).give_proc.recv == 1) {
				(*tp).give_proc.comm_stage = READY;
				(*tp).stage = GIVE_STAGE;
			}
			else {
				(*tp).stage = FINISH;
			}
		}
		break;
	}
	case GIVE_STAGE: {
		MPI_Isend(recv_num, 1, MPI_INT, (*tp).give_proc.rank, num_tag, comm, &(*tp).give_proc.req);
		(*tp).stage = FINISH;
		break;
	}
	default:
		break;
	}
}