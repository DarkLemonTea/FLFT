#ifndef HEADER_FILE
#define HEADER_FILE
#include "ring_head.h"
#endif
#include "ring_count.h"

//=============================================================================================================================
//构造重连进程的链表
//=============================================================================================================================
typedef struct reconnect_Node {
	Comm_proc proc;
	struct time_val {
		double sec;
		double usec;
	} start;

	struct reconnect_Node *next;
}Recon_LL;

typedef struct reconnect_List{
	int need_to_recon_l;
	Recon_LL *left_head;

	int need_to_recon_r;
	Recon_LL *right_head;
	
	int recently_tried_proc; 
	struct time_variable {
		double sec;
		double usec;
	} start;
}Recon_procs;

//创建头结点
Recon_LL *creat_recon_list_head() {
	Recon_LL *head = (Recon_LL*)malloc(sizeof(Recon_LL));
	head->proc.rank = -1;
	head->next = NULL;
	return head;
}

void init_recon_procs(Recon_procs *rp) {
	(*rp).left_head = creat_recon_list_head();
	(*rp).right_head = creat_recon_list_head();

	(*rp).need_to_recon_l = 0;
	(*rp).need_to_recon_r = 0;
}

//查找
int recon_list_in(Recon_LL *head, int target) {
	if (head->next == NULL) {
		return 0;
	}
	else {
		Recon_LL *p;
		p = head->next;
		while (p != NULL) {
			if (p->proc.rank == target) {
				return 1;
			}
			p = p->next;
		}
		return 0;
	}
}

//创建进程
void creat_new_recon_proc(Recon_LL *head, int rank, int stage) {
	Recon_LL *node;
	node = (Recon_LL*)malloc(sizeof(Recon_LL));
	node->proc.rank = rank;
	node->proc.send = stage;
	node->proc.comm_stage = CONNECT_STATE1;

	node->next = head->next;
	head->next = node;
}

//链表删除
void clear_recon_list(Recon_LL *head) {
	Recon_LL *p, *tmp;
	p = head->next;

	while (p != NULL){
		tmp = p;
		p = p->next;
		free(tmp);
	}
	free(head);
}

//寻找最近有效的通信对象
int find_vaild_neighbor(
	int my_rank,
	int comm_size,
	Recon_LL *head,
	int *neighbor_rank
) {
	if (head->next == NULL) { return 0; }

	int res = 0;
	int dis = comm_size;
	Recon_LL *p;
	p = head->next;
	while (p != NULL){
		if (p->proc.comm_stage == CONNECT_FINISH &&
			(my_rank - p->proc.rank + comm_size) % comm_size < dis) {
			dis = (my_rank - p->proc.rank + comm_size) % comm_size;
			*neighbor_rank = p->proc.rank;
			printf("rank %d vaild proc is %d\n", my_rank, *neighbor_rank);
			res = 1;
		}
		p = p->next;
	}
	return res;
}

//=============================================================================================================================
//通信模板
//=============================================================================================================================

//被动通信模板
void passive_comm_template(
	MPI_Comm comm,
	Comm_proc *proc
) {
	switch ((*proc).comm_stage)
	{
	case CONNECT_STATE1: {
		MPI_Iprobe((*proc).rank, RING_CONNECT1, comm, &(*proc).flag, &(*proc).status);
		if ((*proc).flag == 1) {
			MPI_Recv(NULL, 0, MPI_BYTE, (*proc).rank, RING_CONNECT1, comm, MPI_STATUSES_IGNORE);
			MPI_Isend(&(*proc).send, 1, MPI_INT, (*proc).rank, RING_CONNECT2, comm, &(*proc).req);
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
		MPI_Iprobe((*proc).rank, RING_CONNECT3, comm, &(*proc).flag, &(*proc).status);
		if ((*proc).flag == 1) {
			MPI_Recv(&(*proc).recv, 1, MPI_INT, (*proc).rank, RING_CONNECT3, comm, MPI_STATUSES_IGNORE);
			if ((*proc).recv == (*proc).send) {
				MPI_Isend(NULL, 0, MPI_BYTE, (*proc).rank, RING_CONNECT4, comm, &(*proc).req);
				(*proc).comm_stage = CONNECT_STATE4;
			}
			else if ((*proc).recv > (*proc).send) {
				(*proc).comm_stage = CONNECT_LAGGING;
			}
		}
		//MPI_Iprobe((*proc).rank, RING_LAGGING, comm, &(*proc).flag, &(*proc).status);
		//if ((*proc).flag == 1) {
		//	MPI_Recv(NULL, 0, MPI_BYTE, (*proc).rank, RING_LAGGING, comm, MPI_STATUSES_IGNORE);
		//	(*proc).comm_stage = CONNECT_LAGGING;
		//}
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

//主动通信模板
void active_comm_template(
	MPI_Comm comm,
	Comm_proc *proc
) {
	switch ((*proc).comm_stage)
	{
	case CONNECT_STATE1: {
		MPI_Isend(NULL, 0, MPI_BYTE, (*proc).rank, RING_CONNECT1, comm, &(*proc).req);
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
		MPI_Iprobe((*proc).rank, RING_CONNECT2, comm, &(*proc).flag, &(*proc).status);
		if ((*proc).flag == 1) {
			MPI_Recv(&(*proc).recv, 1, MPI_INT, (*proc).rank, RING_CONNECT2, comm, MPI_STATUSES_IGNORE);
			if ((*proc).recv == (*proc).send) {
				MPI_Isend(&(*proc).send, 1, MPI_INT, (*proc).rank, RING_CONNECT3, comm, &(*proc).req);
				(*proc).comm_stage = CONNECT_STATE4;
			}
			else if ((*proc).recv > (*proc).send) {
				(*proc).comm_stage = CONNECT_LAGGING;
			}
		}
		//MPI_Iprobe((*proc).rank, RING_LAGGING, comm, &(*proc).flag, &(*proc).status);
		//if ((*proc).flag == 1) {
		//	MPI_Recv(NULL, 0, MPI_BYTE, (*proc).rank, RING_LAGGING, comm, MPI_STATUSES_IGNORE);
		//	(*proc).comm_stage = CONNECT_LAGGING;
		//}
		break;
	}
	case CONNECT_STATE4: {
		MPI_Test(&(*proc).req, &(*proc).flag, MPI_STATUSES_IGNORE);
		if ((*proc).flag == 1) { (*proc).comm_stage = CONNECT_STATE5; }
		break;
	}
	case CONNECT_STATE5: {
		MPI_Iprobe((*proc).rank, RING_CONNECT4, comm, &(*proc).flag, &(*proc).status);
		if ((*proc).flag == 1) {
			MPI_Recv(NULL, 0, MPI_BYTE, (*proc).rank, RING_CONNECT4, comm, MPI_STATUSES_IGNORE);
			(*proc).comm_stage = CONNECT_FINISH;
		}
		break;
	}
	default:
		break;
	}
}

//被动重连通信模板
void ring_passive_reconnect_template(
	MPI_Comm comm,
	Comm_proc *proc
) {
	switch ((*proc).comm_stage)
	{
	case CONNECT_STATE1: {
		MPI_Isend(&(*proc).send, 1, MPI_INT, (*proc).rank, RING_RECONNECT2, comm, &(*proc).req);
		(*proc).comm_stage = CONNECT_STATE2;
		break;
	}
	case CONNECT_STATE2: {
		MPI_Iprobe((*proc).rank, RING_RECONNECT3, comm, &(*proc).flag, &(*proc).status);
		if ((*proc).flag == 1) {
			MPI_Recv(&(*proc).recv, 1, MPI_INT, (*proc).rank, RING_RECONNECT3, comm, MPI_STATUSES_IGNORE);
			if ((*proc).recv == (*proc).send) {
				MPI_Isend(NULL, 0, MPI_BYTE, (*proc).rank, RING_RECONNECT4, comm, &(*proc).req);
				(*proc).comm_stage = CONNECT_STATE3;
			}
			else if ((*proc).recv > (*proc).send) {
				(*proc).comm_stage = CONNECT_LAGGING;
			}
		}
		//MPI_Iprobe((*proc).rank, RING_LAGGING, comm, &(*proc).flag, &(*proc).status);
		//if ((*proc).flag == 1) {
		//	MPI_Recv(NULL, 0, MPI_BYTE, (*proc).rank, RING_LAGGING, comm, MPI_STATUSES_IGNORE);
		//	(*proc).comm_stage = CONNECT_LAGGING;
		//}
		break;
	}
	case CONNECT_STATE3: {
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

//主动重连通信模板
void ring_active_reconnect_template(
	MPI_Comm comm,
	Comm_proc *proc
) {
	switch ((*proc).comm_stage)
	{
	case CONNECT_STATE1: {
		MPI_Isend(NULL, 0, MPI_BYTE, (*proc).rank, RING_RECONNECT1, comm, &(*proc).req);
		(*proc).comm_stage = CONNECT_STATE2;
		break;
	}
	case CONNECT_STATE2: {
		MPI_Iprobe((*proc).rank, RING_RECONNECT2, comm, &(*proc).flag, &(*proc).status);
		if ((*proc).flag == 1) {
			MPI_Recv(&(*proc).recv, 1, MPI_INT, (*proc).rank, RING_RECONNECT2, comm, MPI_STATUSES_IGNORE);
			if ((*proc).recv == (*proc).send) {
				MPI_Isend(&(*proc).send, 1, MPI_INT, (*proc).rank, RING_RECONNECT3, comm, &(*proc).req);
				(*proc).comm_stage = CONNECT_STATE3;
			}
			else if ((*proc).recv > (*proc).send) {
				(*proc).comm_stage = CONNECT_LAGGING;
			}
		}
		//MPI_Iprobe((*proc).rank, RING_LAGGING, comm, &(*proc).flag, &(*proc).status);
		//if ((*proc).flag == 1) {
		//	MPI_Recv(NULL, 0, MPI_BYTE, (*proc).rank, RING_LAGGING, comm, MPI_STATUSES_IGNORE);
		//	(*proc).comm_stage = CONNECT_LAGGING;
		//}
		break;
	}
	case CONNECT_STATE3: {
		MPI_Iprobe((*proc).rank, RING_RECONNECT4, comm, &(*proc).flag, &(*proc).status);
		if ((*proc).flag == 1) {
			MPI_Recv(NULL, 0, MPI_BYTE, (*proc).rank, RING_RECONNECT4, comm, MPI_STATUSES_IGNORE);
			(*proc).comm_stage = CONNECT_FINISH;
		}
		break;
	}
	default:
		break;
	}
}

//左进程的尝试重连
void ring_passive_reconnect(
	int dec_stage,
	MPI_Comm comm,
	Recon_LL *head
) {
	//探测重连消息
	int flag;
	MPI_Status sta;
	MPI_Iprobe(MPI_ANY_SOURCE, RING_RECONNECT1, comm, &flag, &sta);
	if (flag == 1) {
		if (recon_list_in(head, sta.MPI_SOURCE) == 0) {
			creat_new_recon_proc(head, sta.MPI_SOURCE, dec_stage);
		}
		MPI_Recv(NULL, 0, MPI_BYTE, sta.MPI_SOURCE, RING_RECONNECT1, comm, MPI_STATUSES_IGNORE);
	}

	Recon_LL *p;
	p = head->next;
	while (p != NULL){
		ring_passive_reconnect_template(comm, &(p->proc));
		p = p->next;
	}
}

//右进程的尝试重连
void ring_active_reconnect(
	MPI_Comm comm,
	Recon_LL *head
) {
	Recon_LL *p;
	p = head->next;
	while (p != NULL) {
		ring_active_reconnect_template(comm, &(p->proc));
		p = p->next;
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
void counting_demand(
	int my_rank,
	int comm_size,
	Ring r,
	Recon_procs rp,
	Counter *c
) {
	//相邻进程都成功连接
	if (r.left_proc.comm_stage == CONNECT_FINISH && r.right_proc.comm_stage == CONNECT_FINISH) {
		(*c).demand = 1;
		(*c).stage = STAGE1;
		(*c).r.left_proc.rank = r.left_proc.rank;
		(*c).r.right_proc.rank = r.right_proc.rank;
	}
	else if (r.left_proc.comm_stage != CONNECT_FINISH && r.right_proc.comm_stage == CONNECT_FINISH) {
		(*c).demand = find_vaild_neighbor(my_rank, comm_size, rp.left_head, &((*c).r.left_proc.rank));
		if((*c).demand == 1){
			(*c).r.right_proc.rank = r.right_proc.rank;
			(*c).stage = STAGE1; 
		}
	}
	else if (r.left_proc.comm_stage == CONNECT_FINISH && r.right_proc.comm_stage != CONNECT_FINISH) {
		(*c).demand = find_vaild_neighbor(my_rank, comm_size, rp.right_head, &((*c).r.right_proc.rank));
		if ((*c).demand == 1) {
			(*c).r.left_proc.rank = r.left_proc.rank;
			(*c).stage = STAGE1;
		}
	}
	else {
		(*c).demand = find_vaild_neighbor(my_rank, comm_size, rp.left_head, &((*c).r.left_proc.rank)) &&
			find_vaild_neighbor(my_rank, comm_size, rp.right_head, &((*c).r.right_proc.rank));
		if ((*c).demand == 1) {
			(*c).stage = STAGE1;
		}
	}
}

//判断是否迟到
int is_lagging(
	int my_rank,
	Ring r,
	Recon_procs rp
) {
	int res = 0;

	if (r.left_proc.comm_stage == CONNECT_LAGGING) {
		printf("l proc rank %d tell %d lagging\n", r.left_proc.rank, my_rank);
		res = 1;
		return res;
	}
	if (r.right_proc.comm_stage == CONNECT_LAGGING) { 
		printf("r proc rank %d tell %d lagging\n", r.right_proc.rank, my_rank);
		res = 1; 
		return res;
	}

	Recon_LL *p;
	p = rp.left_head->next;
	while (p != NULL){
		if (p->proc.comm_stage == CONNECT_LAGGING) { 
			printf("l rep rank %d tell %d lagging\n", p->proc.rank, my_rank);
			res = 1; 
			return res;
		}
		p = p->next;
	}

	p = rp.right_head->next;
	while (p != NULL) {
		if (p->proc.comm_stage == CONNECT_LAGGING) {
			printf("r rep rank %d tell %d lagging\n", p->proc.rank, my_rank);
			res = 1;
			return res;
		}
		p = p->next;
	}
	return res;
}

int ring_procs_connect(
	int detector_stage,
	FD_var fd_var,
	MPI_Comm comm,
	P_set lagging_set,
	Ring *ring //用于存储当前有效连接对象
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
	Counter cou;
	init_counter(&cou);
	Comm_proc ver_proc;

	int i = 0, result = -1;
	int T_wait_out = 0;
	struct timeval start, end;
	double cost_time;
	gettimeofday(&start, NULL);
	rp.start.sec = start.tv_sec;
	rp.start.usec = start.tv_usec;
	rp.recently_tried_proc = (*ring).right_proc.rank;

	while (1)
	{
		//环连接
		passive_comm_template(comm, &(*ring).left_proc);
		active_comm_template(comm, &(*ring).right_proc);

		//判断是否重连
		if (T_wait_out == 0) {
			gettimeofday(&end, NULL);
			cost_time = 1000.0 * (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000.0;
			if (cost_time > fd_var.T_wait) {
				T_wait_out = 1;
				if ((*ring).left_proc.comm_stage == CONNECT_FINISH) { rp.need_to_recon_l = 0; }
				else{ rp.need_to_recon_l = 1; }

				if ((*ring).right_proc.comm_stage == CONNECT_FINISH) { rp.need_to_recon_r = 0; }
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
					rp.recently_tried_proc = find_valid(comm_size, rp.recently_tried_proc,
						lagging_set.num, lagging_set.procs, 'r');
					//将该进程加入链表中
					creat_new_recon_proc(rp.right_head, rp.recently_tried_proc, detector_stage);
				}
				ring_active_reconnect(comm, rp.right_head);
			}
			//如果直连进程到达，就不需要再重连，关闭重连通道
			if ((*ring).left_proc.comm_stage == CONNECT_FINISH) { rp.need_to_recon_l = 0; }
			if ((*ring).right_proc.comm_stage == CONNECT_FINISH) { rp.need_to_recon_r = 0; }
		}

		//判断参与计数的进程
		if (cou.demand == 0) { counting_demand(my_rank, comm_size, *ring, rp, &cou); }
		//有计数需求，执行计数
		if (cou.demand == 1) { arrived_procs_count(comm, my_rank, &cou); }
		//如果计数完成
		if (cou.stage == FINISH) {
			//printf("rank %d count num is %d\n", my_rank, cou.sum);

			if (cou.sum < min_arr_num) {
				//计数未达标，初始化counter，重新计数
				init_counter(&cou);
			}
			else {
				gettimeofday(&end, NULL);
				cost_time = 1000.0 * (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000.0;

				//达到计数条件，通过
				(*ring).left_proc.rank = cou.r.left_proc.rank;
				(*ring).right_proc.rank = cou.r.right_proc.rank;
				if (cou.sum == comm_size - lagging_set.num) {
					result = RING_CONNECT_SUCCESS;
				}
				else {
					result = RING_CONNECT_THROUGH;
				}		
				break;
			}
		}
		
		//判断是否迟到
		if (is_lagging(my_rank,*ring, rp)) {
			result = RING_CONNECT_LAGGING;
		}

		//判断是否有错误重连
		/*MPI_Iprobe(MPI_ANY_SOURCE, RING_MISCONNECTION, comm, &ver_proc.flag, &ver_proc.status);
		if (ver_proc.flag == 1) {
			MPI_Recv(NULL, 0, MPI_BYTE, ver_proc.status.MPI_SOURCE, RING_MISCONNECTION, comm, MPI_STATUSES_IGNORE);
			result = RING_CONNECT_LAGGING;
		}*/
		//如果迟到，且右边进程有正确连接，发送通知
		if (result == RING_CONNECT_LAGGING) {
			/*if (cou.r.right_proc.rank > -1) {
				MPI_Isend(NULL, 0, MPI_BYTE, cou.r.right_proc.rank, RING_MISCONNECTION, comm, &cou.r.right_proc.req);
			}*/
			break;
		}

		gettimeofday(&end, NULL);
		cost_time = 1000.0 * (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000.0;
		//如果超过最大时限，返回失败
		if (cost_time > fd_var.T_max) {
			result = RING_CONNECT_FAILURE;
			break;
		}
	}

	//删除链表
	clear_recon_list(rp.left_head);
	clear_recon_list(rp.right_head);

	return result;
}


//=============================================================================================================================
//lagging进程的主动阻塞状态，被捡回后需要确定一个rescuer
//=============================================================================================================================
void block_and_rescuer_confirm(
	MPI_Comm comm,
	FD_var fd,
	int lagging_num,
	int *lagging_procs,
	Comm_proc *rescuer
) {
	int my_rank, comm_size;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &comm_size);
	printf("rank %d block!\n", my_rank);

	double block_time = 5 * fd.T_wait;
	Comm_proc savior, tmp_proc;
	savior.comm_stage = CONNECT_LAGGING;
	double cost_time;
	struct timeval st,start,time;
	
	gettimeofday(&st, NULL);
	
	while (savior.comm_stage != CONNECT_LAGGING_DONE) {
		switch (savior.comm_stage)
		{
		case CONNECT_LAGGING: {
			//进程响应后，应该将进程加入待捡回进程
			savior.comm_stage = CONNECT_LAGGING_STATE1;
			gettimeofday(&start, NULL);
			break;
		}
		case CONNECT_LAGGING_STATE1: {
			MPI_Iprobe(MPI_ANY_SOURCE, RING_LAGGING_PROBE, comm, &savior.flag, &savior.status);
			if (savior.flag == 1) {
				printf("recv lagging probe from rank %d\n", savior.status.MPI_SOURCE);
				savior.rank = savior.status.MPI_SOURCE;
				MPI_Recv(NULL, 0, MPI_BYTE, savior.rank, RING_LAGGING_PROBE, comm, MPI_STATUSES_IGNORE);
				MPI_Isend(NULL, 0, MPI_BYTE, savior.rank, RING_LAGGING_RESPOND, comm, &savior.req);
				savior.comm_stage = CONNECT_LAGGING_STATE2;
			}
			else {
				gettimeofday(&time, NULL);
				cost_time = 1000.0 * (time.tv_sec - start.tv_sec) + (time.tv_usec - start.tv_usec) / 1000.0;
				//如果超过阻塞时间，寻找下一个有效的对象求助
				if (cost_time > block_time) {
					savior.comm_stage = CONNECT_LAGGING_TIMEOUT;		
				}
			}
			break;
		}
		case CONNECT_LAGGING_STATE2: {
			MPI_Test(&savior.req, &savior.flag, MPI_STATUSES_IGNORE);
			if (savior.flag == 1) {
				savior.comm_stage = CONNECT_LAGGING_STATE3;
			}
			break;
		}
		case CONNECT_LAGGING_STATE3: {
			MPI_Iprobe(savior.rank, RING_LAGGING_CONFIRM, comm, &savior.flag, MPI_STATUSES_IGNORE);
			if (savior.flag == 1) {
				(*rescuer).rank = savior.rank;
				(*rescuer).comm_stage = READY;
				MPI_Recv(&(*rescuer).recv, 1, MPI_INT, (*rescuer).rank, RING_LAGGING_CONFIRM, comm, MPI_STATUSES_IGNORE);
				savior.comm_stage = CONNECT_LAGGING_DONE;	
				printf("rank %d be rescued by %d\n", my_rank, (*rescuer).rank);
			}
			else
			{
				gettimeofday(&time, NULL);
				cost_time = 1000.0 * (time.tv_sec - start.tv_sec) + (time.tv_usec - start.tv_usec) / 1000.0;
				//如果超过阻塞时间，寻找下一个有效的对象求助
				if (cost_time > block_time) {
					savior.comm_stage = CONNECT_LAGGING_TIMEOUT;
				}
			}
			break;
		}
		case CONNECT_LAGGING_TIMEOUT: {
			//随机寻找一个有效进程发送消息
			tmp_proc.rank = get_rand_proc(my_rank, comm_size);
			tmp_proc.comm_stage = CONNECT_STATE1;
			gettimeofday(&start, NULL);
			while (1) {
				ring_active_reconnect_template(comm, &tmp_proc);
				if (tmp_proc.comm_stage == CONNECT_LAGGING) {
					break;
				}
				else {
					gettimeofday(&time, NULL);
					cost_time = 1000.0 * (time.tv_sec - start.tv_sec) + (time.tv_usec - start.tv_usec) / 1000.0;
					if (cost_time > block_time) {
						tmp_proc.rank = get_rand_proc(my_rank, comm_size);
						tmp_proc.comm_stage = CONNECT_STATE1;
						gettimeofday(&start, NULL);
					}
				}
			}

			savior.rank = tmp_proc.rank;
			savior.comm_stage = CONNECT_LAGGING_STATE1;
			gettimeofday(&start, NULL);
			break;
		}
		default:
			break;
		}

		gettimeofday(&time, NULL);
		cost_time = 1000.0 * (time.tv_sec - st.tv_sec) + (time.tv_usec - st.tv_usec) / 1000.0;
		//如果超过阻塞时间，寻找下一个有效的对象求助
		if (cost_time > fd.T_max) {
			(*rescuer).comm_stage = CONNECT_TIMEOUT_BREAK;
			break;
		}
	}
}