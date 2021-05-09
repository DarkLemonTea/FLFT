#pragma once
#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>
#include<math.h>

#include<sys/time.h>
#include<unistd.h>

#define SUCCESS 0
#define ERROR 10
#define EMPTY -1

//=============================================================================================================================
//标识管理
//=============================================================================================================================

#define READY 0
#define STAGE1 1
#define STAGE2 2
#define STAGE3 3
#define STAGE4 4
#define FINISH 5

#define RECV_STAGE 6
#define REQ_CONFIRM_STAGE 7
#define REQ_RECV_STAGE 8
#define GIVE_STAGE 9

//本地操作执行阶段
#define CONNECT_STATE1 1
#define CONNECT_STATE2 2
#define CONNECT_STATE3 3
#define CONNECT_STATE4 4
#define CONNECT_STATE5 5
#define CONNECT_FINISH 6
#define CONNECT_LAGGING 7
#define CONNECT_LAGGING_STATE1 8
#define CONNECT_LAGGING_STATE2 9
#define CONNECT_LAGGING_STATE3 10
#define CONNECT_LAGGING_RETRIEVE 11
#define CONNECT_LAGGING_DONE 13
#define CONNECT_LAGGING_TIMEOUT 14
#define CONNECT_TIMEOUT_BREAK 15
#define CONNECT_FAILURE 16

//=============================================================================================================================
//基础结构
//=============================================================================================================================

typedef struct comm_proc {
	int rank;
	int flag;
	int comm_stage; //连接阶段
	MPI_Request req;
	MPI_Status status;
	int send;
	int recv;
}Comm_proc;

void init_proc(
	int rank,
	Comm_proc *proc
) {
	proc->rank = rank;
	proc->comm_stage = READY;
}

//=============================================================================================================================
//链表，用于进程的捡回操作
//=============================================================================================================================

typedef struct revive_Node {
	Comm_proc revive_proc;
	struct revive_Node* next;
}Revive_LL;

//用于存放集合的指针
typedef struct Procs_set {
	int num;
	int *procs;
}P_set;

void init_p_set(P_set *p) {
	(*p).num = 0;
	(*p).procs = NULL;
}

//=============================================================================================================================
//基础操作
//=============================================================================================================================

//查找是否在列表中
int in(P_set l, int target) {
	if (l.num == 0) return 0;
	else {
		int i;
		for (i = 0; i < l.num; i++) {
			if (l.procs[i] == target) return 1;
		}
		return 0;
	}
}

//合并数组
void merge_procs(
	//input
	P_set a,
	P_set b,
	//output
	P_set *c
) {
	int num = 0, i;
	int *procs;
	if (a.num == 0 && b.num == 0) {
		(*c).num = num;
		(*c).procs = NULL;
	}
	else if (b.num == 0) {
		num = a.num;
		procs = a.procs;
	}
	else if (a.num == 0) {
		num = b.num;
		procs = b.procs;
	}
	else {
		num = a.num + b.num;
		procs = (int*)calloc(num, sizeof(int));
		for (i = 0; i < num; i++) {
			if (i < b.num) {
				procs[i] = b.procs[i];
			}
			else {
				procs[i] = a.procs[i - b.num];
			}
		}
	}
	(*c).num = num;
	(*c).procs = procs;
}

//合并三数组
void merge_tri_procs(
	//input
	P_set a,
	P_set b,
	P_set c,
	//output
	P_set *d
) {
	P_set t;
	merge_procs(a,b,&t);
	merge_procs(c,t,d);
}

//随机选一个进程
int get_rand_proc(
	int my_rank,
	int comm_size
) {
	int proc = rand() % comm_size;
	while (proc == my_rank ||
		proc == (my_rank + comm_size - 1) % comm_size ||
		proc == (my_rank + 1) % comm_size) {
		proc = rand() % comm_size;
	}
	return proc;
}


//=============================================================================================================================
//数据发送模板
//=============================================================================================================================

//发送数组模板
void send_procs_template(
	MPI_Comm comm,
	Comm_proc *proc,
	int num_tag,
	int procs_tag,
	P_set s
) {
	switch ((*proc).comm_stage)
	{
	case READY:
		MPI_Isend(&(s.num), 1, MPI_INT, (*proc).rank, num_tag, comm, &(*proc).req);
		(*proc).comm_stage = STAGE1;
		break;
	case STAGE1:
		MPI_Test(&(*proc).req, &(*proc).flag, MPI_STATUSES_IGNORE);
		if ((*proc).flag == 1) {
			if (s.num > 0) {
				MPI_Isend(s.procs, s.num, MPI_INT, (*proc).rank, procs_tag, comm, &(*proc).req);
				//printf("send[0] is %d\n", send_procs[0]);
				(*proc).comm_stage = STAGE2;
			}
			else {
				(*proc).comm_stage = FINISH;
			}
		}
		break;
	case STAGE2:
		MPI_Test(&(*proc).req, &(*proc).flag, MPI_STATUSES_IGNORE);
		if ((*proc).flag == 1) {
			(*proc).comm_stage = FINISH;
		}
		break;
	default:
		break;
	}
}

//接收数组模板
void recv_procs_template(
	MPI_Comm comm,
	Comm_proc *proc,
	int num_tag,
	int procs_tag,
	P_set *r
) {
	switch ((*proc).comm_stage)
	{
	case READY: {
		MPI_Iprobe((*proc).rank, num_tag, comm, &(*proc).flag, MPI_STATUSES_IGNORE);
		if ((*proc).flag == 1) {
			MPI_Recv(&((*r).num), 1, MPI_INT, (*proc).rank, num_tag, comm, MPI_STATUSES_IGNORE);
			if ((*r).num > 0) {
				(*r).procs = (int*)calloc((*r).num, sizeof(int));
				MPI_Irecv((*r).procs, (*r).num, MPI_INT, (*proc).rank, procs_tag, comm, &(*proc).req);
				(*proc).comm_stage = STAGE1;
			}
			else {
				(*proc).comm_stage = FINISH;
			}
		}
		break;
	}
	case STAGE1: {
		MPI_Test(&(*proc).req, &(*proc).flag, MPI_STATUSES_IGNORE);
		if ((*proc).flag == 1) {
			(*proc).comm_stage = FINISH;
		}
		break;
	}
	default:
		break;
	}
}

void send_num_template(
	MPI_Comm comm,
	Comm_proc *proc,
	int num_tag,
	int send_num
) {
	switch ((*proc).comm_stage)
	{
	case READY:
		if ((*proc).rank == EMPTY) {
			(*proc).comm_stage = FINISH;
			break;
		}
		MPI_Isend(&send_num, 1, MPI_INT, (*proc).rank, num_tag, comm, &(*proc).req);
		(*proc).comm_stage = STAGE1;
		break;
	case STAGE1:
		MPI_Test(&(*proc).req, &(*proc).flag, MPI_STATUSES_IGNORE);
		if ((*proc).flag == 1) {
			(*proc).comm_stage = FINISH;
		}
		break;
	default:
		break;
	}
}

//=============================================================================================================================
//排序
//=============================================================================================================================
void swap(int arr[], int low, int high)
{
	int temp;
	temp = arr[low];
	arr[low] = arr[high];
	arr[high] = temp;
}

int Partition(int array[], int low, int high) {
	int base = array[low];
	while (low < high) {
		while (low < high && array[high] >= base) {
			high--;
		}
		swap(array, low, high);//array[low] = array[high];
		while (low < high && array[low] <= base) {
			low++;
		}
		swap(array, low, high);//array[high] = array[low];
	}
	array[low] = base;
	return low;
}

void QuickSort(int array[], int low, int high) {
	if (low < high) {
		int base = Partition(array, low, high);
		QuickSort(array, low, base - 1);
		QuickSort(array, base + 1, high);
	}
}

