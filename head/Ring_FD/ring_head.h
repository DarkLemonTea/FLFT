#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>
#include<math.h>
#include<sys/time.h>
#include<unistd.h>

//=============================================================================================================================
//标识管理
//=============================================================================================================================

//————————————————————————————————————————————————————————
//返回结果
#define RING_CONNECT_SUCCESS 0
#define RING_CONNECT_THROUGH 1
#define RING_CONNECT_FAILURE 2
#define RING_CONNECT_LAGGING 3

#define FD_SUCCESS 0
#define FD_FAILURE 1
#define FD_REVIVE 2
#define FD_PASS 3

#define READY 0
#define STAGE1 1
#define STAGE2 2
#define STAGE3 3
#define STAGE4 4
#define FINISH 5

#define GATHER_SUCCESS 0
#define CONFIRM_LAGGING 1
#define UNKNOWN_ERROR 2

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

#define EMPTY -1

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

//————————————————————————————————————————————————————————
//容错集合通信
#define TREE_BCAST 30005
#define BRUCK_ALLGATHER 30006
#define BRUCK_RECVNUM 30007
#define BRUCK_MULTIGATHER 30008


//=============================================================================================================================
//基础结构和操作
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

//查找是否在列表中
int in(int num, int *list, int target) {
	if (num == 0) return 0;
	else {
		int i;
		for (i = 0; i < num; i++) {
			if (list[i] == target) return 1;
		}
		return 0;
	}
}

//合并数组
void merge_procs(
	//input
	int a_num,
	int b_num,
	int *a_procs,
	int *b_procs,
	//output
	int *c_num,
	int **c_procs
) {
	int num = 0, i;
	int *procs;
	if (a_num == 0 && b_num == 0) {
		*c_num = num;
		*c_procs = NULL;
	}
	else if (b_num == 0) {
		num = a_num;
		procs = a_procs;
	}
	else if (a_num == 0) {
		num = b_num;
		procs = b_procs;
	}
	else {
		num = a_num + b_num;
		procs = (int*)calloc(num, sizeof(int));
		for (i = 0; i < num; i++) {
			if (i < b_num) {
				procs[i] = b_procs[i];
			}
			else {
				procs[i] = a_procs[i - b_num];
			}
		}
	}
	*c_num = num;
	*c_procs = procs;
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
//Ring
//=============================================================================================================================

typedef struct ring {
	Comm_proc left_proc;
	Comm_proc right_proc;
}Ring;

int left_proc(int comm_size, int my_rank) { return (comm_size + my_rank - 1) % comm_size; }
int right_proc(int comm_size, int my_rank) { return (my_rank + 1) % comm_size; }

//寻找最近有效对象
int find_valid(int comm_size, int rank, int lagging_num, int *lagging_procs, char l_or_r) {
	int tmp_proc;
	if (l_or_r == 'l') {
		tmp_proc = left_proc(comm_size, rank);
		while (in(lagging_num, lagging_procs, tmp_proc)) { tmp_proc = left_proc(comm_size, tmp_proc); }		
	}
	else if (l_or_r == 'r') {
		tmp_proc = right_proc(comm_size, rank);
		while (in(lagging_num, lagging_procs, tmp_proc)) { tmp_proc = right_proc(comm_size, tmp_proc); }	
	}
	return tmp_proc;
}

//初始化环
void init_ring(
	int comm_size,
	int my_rank,
	int detector_stage,
	int lagging_num,
	int *lagging_procs,
	Ring *ring
) {
	ring->left_proc.rank = find_valid(comm_size, my_rank, lagging_num, lagging_procs, 'l');
	ring->right_proc.rank = find_valid(comm_size, my_rank, lagging_num, lagging_procs, 'r');
	ring->left_proc.comm_stage = CONNECT_STATE1;
	ring->right_proc.comm_stage = CONNECT_STATE1;
	ring->left_proc.send = detector_stage;
	ring->right_proc.send = detector_stage;
}


//发送数组模板
void send_procs_template(
	MPI_Comm comm,
	Comm_proc *proc,
	int num_tag,
	int procs_tag,
	int send_num,
	int *send_procs
) {
	switch ((*proc).comm_stage)
	{
	case READY:
		MPI_Isend(&send_num, 1, MPI_INT, (*proc).rank, num_tag, comm, &(*proc).req);
		(*proc).comm_stage = STAGE1;
		break;
	case STAGE1:
		MPI_Test(&(*proc).req, &(*proc).flag, MPI_STATUSES_IGNORE);
		if ((*proc).flag == 1) {
			if (send_num > 0) {
				MPI_Isend(send_procs, send_num, MPI_INT, (*proc).rank, procs_tag, comm, &(*proc).req);
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
	int *recv_num,
	int **recv_procs
) {
	switch ((*proc).comm_stage)
	{
	case READY: {
		MPI_Iprobe((*proc).rank, num_tag, comm, &(*proc).flag, MPI_STATUSES_IGNORE);
		if ((*proc).flag == 1) {
			MPI_Recv(recv_num, 1, MPI_INT, (*proc).rank, num_tag, comm, MPI_STATUSES_IGNORE);
			if (*recv_num > 0) {
				*recv_procs = (int*)calloc(*recv_num, sizeof(int));
				MPI_Irecv(*recv_procs, *recv_num, MPI_INT, (*proc).rank, procs_tag, comm, &(*proc).req);
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
			//printf("recv procs[0] is %d\n", (*recv_procs)[0]);
			(*proc).comm_stage = FINISH;
		}
		break;
	}
	default:
		break;
	}
}

void statistics_lagging_procs(
	//input
	Ring ring,
	int my_rank,
	int comm_size,
	//output
	int *local_num,
	int **local_procs
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
		*local_num = 0;
		*local_procs = NULL;
	}
	else {
		*local_num = num;
		procs = (int*)calloc(num, sizeof(int));
		for (i = 0; i < num; i++) {
			procs[i] = (ring.left_proc.rank + 1 + i) % comm_size;
		}
		*local_procs = procs;
	}
}

//===========================================================================================================================

//获得最后的进程号
int get_the_last(
	int comm_size,
	int recv_num,
	int *recv_procs
) {
	int last_proc = comm_size - 1;
	while (in(recv_num, recv_procs, last_proc)) {
		last_proc -= 1;
	}
	return last_proc;
}

//根据进程号和L，获取坐标，用于通信拓扑
int get_comm_index(
	int my_rank,
	int lagging_num,
	int *lagging_procs,
	int comm_size
) {
	int i, ind = my_rank;
	for (i = 0; i < lagging_num; i++) {
		if (lagging_procs[i] < my_rank) {
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
	int lagging_num,
	int *lagging_procs,
	int comm_size
) {
	//已排序
	int i = 0;
	int real_num = index;
	while (i < lagging_num)
	{
		if (lagging_procs[i] <= real_num) {
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

//快排
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

typedef struct Set_pointers {
	int current_stage;
	P_set Lagging_procs;         //迟到进程列表
	Revive_LL *local_re_phead;   //存储恢复进程的临时链表
	P_set Revive_procs;          //恢复进程列表
	Ring ring;                   //环结构
}Detector;

void init_detector(Detector *sp) {
	//初始化
	(*sp).Lagging_procs.num = 0;
	(*sp).Lagging_procs.procs = NULL;
	
	//创建一个头结点，不代表任何数据
	(*sp).local_re_phead = (Revive_LL*)malloc(sizeof(Revive_LL));
	(*sp).local_re_phead->revive_proc.rank = -1;
	(*sp).local_re_phead->next = NULL;

	(*sp).Revive_procs.num = 0;
	(*sp).Revive_procs.procs = NULL;
}

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
		while (p != NULL){
			if (p->revive_proc.rank == target) {
				printf("rank %d in lagging procs\n", p->revive_proc.rank);
				return 1;
			}
			p = p->next;
		}
		return 0;
	}
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
