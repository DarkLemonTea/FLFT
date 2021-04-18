#ifndef HEADER_FILE
#define HEADER_FILE
#include"ring_head.h"
#include"ring_connect.h"
#endif

//=============================================================================================================================
//(1)探测消息并捡回
//=============================================================================================================================

void probe_connect1(
	MPI_Comm comm,
	int probe_tag,
	int respond_tag,
	int stage,
	Detector *sp
) {
	Comm_proc prober;
	prober.flag = 1;
	while (prober.flag == 1){
		MPI_Iprobe(MPI_ANY_SOURCE, probe_tag, comm, &prober.flag, &prober.status);

		if (prober.flag == 1) {
			if (in((*sp).Lagging_procs.num, (*sp).Lagging_procs.procs, prober.status.MPI_SOURCE) &&
				!linked_list_in((*sp).local_re_phead, prober.status.MPI_SOURCE)) {

				//printf("probe message from %d\n", prober.status.MPI_SOURCE);
				creat_and_init((*sp).local_re_phead, prober.status.MPI_SOURCE, stage);

				MPI_Recv(NULL, 0, MPI_BYTE, (*sp).local_re_phead->next->revive_proc.rank,
					probe_tag, comm, MPI_STATUSES_IGNORE);
				MPI_Isend(&((*sp).local_re_phead->next->revive_proc.send), 1, MPI_INT,
					(*sp).local_re_phead->next->revive_proc.rank, respond_tag, comm,
					&((*sp).local_re_phead->next->revive_proc.req));

				//(*sp).loc_re_p.head->next->revive_proc.send = stage;
				//(*sp).loc_re_p.head->next->revive_proc.comm_stage = STAGE1;
			}
			else if (linked_list_in((*sp).local_re_phead, prober.status.MPI_SOURCE)) {
				MPI_Recv(NULL, 0, MPI_BYTE, prober.status.MPI_SOURCE, probe_tag, comm, MPI_STATUSES_IGNORE);
			}
			else {
				break;
			}
		}
	}
}

void probe_connect2(
	MPI_Comm comm,
	int probe_tag,
	int respond_tag,
	int stage,
	Detector *sp
) {
	int buff;
	Comm_proc prober;
	prober.flag = 1;
	while (prober.flag == 1) {
		MPI_Iprobe(MPI_ANY_SOURCE, probe_tag, comm, &prober.flag, &prober.status);
		if (prober.flag == 1) {
			if (in((*sp).Lagging_procs.num, (*sp).Lagging_procs.procs, prober.status.MPI_SOURCE) &&
				!linked_list_in((*sp).local_re_phead, prober.status.MPI_SOURCE)) {

				//printf("probe message from %d\n", prober.status.MPI_SOURCE);
				creat_and_init((*sp).local_re_phead, prober.status.MPI_SOURCE, stage);

				MPI_Recv(&prober.recv, 1, MPI_INT, (*sp).local_re_phead->next->revive_proc.rank,
					probe_tag, comm, MPI_STATUSES_IGNORE);
				MPI_Isend(&((*sp).local_re_phead->next->revive_proc.send), 1, MPI_INT,
					(*sp).local_re_phead->next->revive_proc.rank, respond_tag, comm,
					&((*sp).local_re_phead->next->revive_proc.req));

				//(*sp).loc_re_p.head->next->revive_proc.send = stage;
				//(*sp).loc_re_p.head->next->revive_proc.comm_stage = STAGE1;
			}
			else if (linked_list_in((*sp).local_re_phead, prober.status.MPI_SOURCE)) {
				MPI_Recv(&buff, 1, MPI_INT, prober.status.MPI_SOURCE, probe_tag, comm, MPI_STATUSES_IGNORE);
			}
		}
	}
}

//检测函数，如果检测到lagging进程消息，回复，并添加到链表中
int probe_revive_procs(
	int detector_stage,
	MPI_Comm comm,
	Detector *sp
) {
	if ((*sp).Lagging_procs.num == 0) {
		return 0;
	}

	/*int my_rank;
	MPI_Comm_rank(comm, &my_rank);
	printf("rank %d lagging num is %d\n", my_rank, (*sp).Lagging_procs.num);*/

	//Tag为RING_CONNECT1的消息进程
	probe_connect1(comm, RING_CONNECT1, RING_CONNECT2, detector_stage, sp);

	//Tag为RING_CONNECT2的消息进程
	probe_connect2(comm, RING_CONNECT2, RING_CONNECT3, detector_stage, sp);

	//Tag为RING_RECONNECT1的消息进程
	probe_connect1(comm, RING_RECONNECT1, RING_RECONNECT2, detector_stage, sp);

	//Tag为RING_RECONNECT2的消息进程
	probe_connect2(comm, RING_RECONNECT2, RING_RECONNECT3, detector_stage, sp);

	return 0;
}

//=============================================================================================================================
//(2)确认local进程，与进程通信，同时多聚集捡回的进程
//=============================================================================================================================
void rescue_lagging_procs(
	int stage,
	MPI_Comm comm,
	Comm_proc *proc
) {
	switch ((*proc).comm_stage)
	{
	case STAGE1: {
		MPI_Isend(NULL, 0, MPI_BYTE, (*proc).rank, RING_LAGGING_PROBE, comm, &(*proc).req);
		(*proc).comm_stage = STAGE2;
		break;
	}
	case STAGE2: {
		MPI_Test(&(*proc).req, &(*proc).flag, MPI_STATUSES_IGNORE);
		if ((*proc).flag == 1) {
			(*proc).comm_stage = STAGE3;
		}
		break;
	}
	case STAGE3: {
		MPI_Iprobe((*proc).rank, RING_LAGGING_RESPOND, comm, &(*proc).flag, MPI_STATUSES_IGNORE);
		if ((*proc).flag == 1) {
			//printf("rescuer recv respond from %d\n", (*proc).rank);
			MPI_Recv(NULL, 0, MPI_BYTE, (*proc).rank, RING_LAGGING_RESPOND, comm, MPI_STATUSES_IGNORE);
			MPI_Isend(&stage, 1, MPI_INT, (*proc).rank, RING_LAGGING_CONFIRM, comm, &(*proc).req);
			(*proc).comm_stage = STAGE4;
		}
		break;
	}
	case STAGE4: {
		MPI_Test(&(*proc).req, &(*proc).flag, MPI_STATUSES_IGNORE);
		if ((*proc).flag == 1) {
			(*proc).comm_stage = FINISH;
			printf("rescue rank %d\n", (*proc).rank);
		}
		break;
	}
	case FINISH: {
		break;
	}
	default:
		break;
	}
}

void ring_multigather_revive_procs(
	MPI_Comm comm,
	Ring ring,
	int local_revive_num,
	int *local_revive_procs,
	//output
	int *ring_revive_procs_num,
	int **ring_revive_procs_procs
) {
	int my_rank, comm_size;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &comm_size);

	int ide, i;
	if (ring.left_proc.rank > my_rank) { ide = 0; }//起始进程
	else if (my_rank > ring.right_proc.rank) { ide = 2; }//末尾进程
	else { ide = 1; }

	int recv_num = 0, send_num = 0, total_num = 0;
	int *recv_procs = NULL, *send_procs = NULL, *total_procs = NULL;

	Comm_proc left_proc, right_proc;
	init_proc(ring.left_proc.rank, &left_proc);
	init_proc(ring.right_proc.rank, &right_proc);

	switch (ide)
	{
	case 0:
	{	//第一轮环状多播，非阻塞实现	
		while (1) {
			send_procs_template(comm, &right_proc, RING_GATHER_REVIVE_NUM, RING_GATHER_REVIVE_PROCS, local_revive_num, local_revive_procs);
			if (right_proc.comm_stage == FINISH) { break; }
		}

		//第二轮按照环形
		left_proc.comm_stage = READY; right_proc.comm_stage = READY;
		while (1)
		{
			//接收数组
			if (left_proc.comm_stage != FINISH) {
				recv_procs_template(comm, &left_proc, RING_MULTICAST_REVIVE_NUM, RING_MULTICAST_REVIVE_PROCS, &total_num, &total_procs);
				if (left_proc.comm_stage == FINISH) {
					*ring_revive_procs_num = total_num;
					*ring_revive_procs_procs = total_procs;
				}
			}
			else {
				//发送数组
				send_procs_template(comm, &right_proc, RING_MULTICAST_REVIVE_NUM, RING_MULTICAST_REVIVE_PROCS, total_num, total_procs);
			}
			if (right_proc.comm_stage == FINISH) { break; }
		}
		break;
	}
	case 1:
	{
		//第一轮	
		while (1) {
			if (left_proc.comm_stage != FINISH) {
				recv_procs_template(comm, &left_proc, RING_GATHER_REVIVE_NUM, RING_GATHER_REVIVE_PROCS, &recv_num, &recv_procs);
			}
			else {
				if (right_proc.comm_stage == READY) {
					merge_procs(local_revive_num, recv_num, local_revive_procs, recv_procs, &send_num, &send_procs);
					send_procs_template(comm, &right_proc, RING_GATHER_REVIVE_NUM, RING_GATHER_REVIVE_PROCS, send_num, send_procs);
				}
				else {
					send_procs_template(comm, &right_proc, RING_GATHER_REVIVE_NUM, RING_GATHER_REVIVE_PROCS, send_num, send_procs);
				}
			}
			if (right_proc.comm_stage == FINISH) { break; }
		}

		//第二轮按照环形
		left_proc.comm_stage = READY; right_proc.comm_stage = READY;
		while (1)
		{
			//接收数组
			if (left_proc.comm_stage != FINISH) {
				recv_procs_template(comm, &left_proc, RING_MULTICAST_REVIVE_NUM, RING_MULTICAST_REVIVE_PROCS, &total_num, &total_procs);
				if (left_proc.comm_stage == FINISH) {
					*ring_revive_procs_num = total_num;
					*ring_revive_procs_procs = total_procs;
				}
			}
			else {
				//发送数组
				send_procs_template(comm, &right_proc, RING_MULTICAST_REVIVE_NUM, RING_MULTICAST_REVIVE_PROCS, total_num, total_procs);
			}
			if (right_proc.comm_stage == FINISH) { break; }
		}
		break;
	}
	case 2:
	{
		while (1) {
			recv_procs_template(comm, &left_proc, RING_GATHER_REVIVE_NUM, RING_GATHER_REVIVE_PROCS, &recv_num, &recv_procs);
			if (left_proc.comm_stage == FINISH) { break; }
		}

		merge_procs(local_revive_num, recv_num, local_revive_procs, recv_procs, &total_num, &total_procs);
		*ring_revive_procs_num = total_num;
		*ring_revive_procs_procs = total_procs;

		//第二轮，多播
		left_proc.comm_stage = READY; right_proc.comm_stage = READY;
		while (1)
		{
			send_procs_template(comm, &right_proc, RING_MULTICAST_REVIVE_NUM, RING_MULTICAST_REVIVE_PROCS, total_num, total_procs);
			if (right_proc.comm_stage == FINISH) { break; }
		}
		break;
	}
	default: { break; }
	}
}

void remove_revive_procs(
	Detector *sp
) {
	int num = (*sp).Lagging_procs.num - (*sp).Revive_procs.num;
	
	if (num == 0) {
		(*sp).Lagging_procs.num = 0;
		(*sp).Lagging_procs.procs = NULL;	
	}
	else {
		int *procs;
		procs = (int*)calloc(num, sizeof(int));
		int ind = 0, i;
		for (i = 0; i < (*sp).Lagging_procs.num; i++) {
			if (in((*sp).Revive_procs.num, (*sp).Revive_procs.procs, (*sp).Lagging_procs.procs[i])) {
				continue;
			}
			else {
				procs[i] = (*sp).Lagging_procs.procs[i];
				ind += 1;
			}
		}
		(*sp).Lagging_procs.num = num;
		(*sp).Lagging_procs.procs = procs;
	}
}

//聚集捡回进程
void ring_retrieve_procs(
	int detector_stage,
	MPI_Comm comm,
	double T_retrieve,
	Detector *sp
) {
	int my_rank, comm_size;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &comm_size);
	//init_ring(comm_size, my_rank, detector_stage, (*sp).Lagging_procs.num, (*sp).Lagging_procs.procs, &(*sp).ring);
	//printf("rank %d left is %d,right is %d\n", my_rank, (*sp).ring.left_proc.rank, (*sp).ring.right_proc.rank);

	(*sp).Last_lagging_procs = (*sp).Lagging_procs;

	(*sp).Revive_procs.num = 0;
	(*sp).Revive_procs.procs = NULL;

	//第一步：捡回进程，通过通信确认
	struct timeval start, end;
	double cost_time = 0;
	gettimeofday(&start, NULL);
	Revive_LL *p, *pre;
	while (cost_time < T_retrieve)
	{
		probe_revive_procs(detector_stage, comm, sp);

		p = (*sp).local_re_phead->next;
		while (p != NULL){
			rescue_lagging_procs(detector_stage, comm, &(p->revive_proc));
			p = p->next;
		}
		gettimeofday(&end, NULL);
		cost_time = 1000.0 * (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000.0;
	}

	//第二步：筛选出捡回的本地进程，捡回失败的删除
	int local_revive_num = 0;
	int *local_revive_procs;
	
	//第一次遍历，计数
	pre = (*sp).local_re_phead;
	p = (*sp).local_re_phead->next;
	while (p != NULL) {
		if (p->revive_proc.comm_stage == FINISH) {
			local_revive_num += 1;
			p->revive_proc.comm_stage = READY;
			pre = p;
			p = p->next;
		}
		else {
			pre->next = p->next;
			free(p);
			p = pre->next;
		}
	}
	local_revive_procs = (int*)calloc(local_revive_num, sizeof(int));

	//第二次遍历，填入
	pre = (*sp).local_re_phead;
	p = (*sp).local_re_phead->next;
	int i = 0;
	while (p != NULL) {
		local_revive_procs[i] = p->revive_proc.rank;
		i += 1;
		pre = p;
		p = p->next;
	}

	ring_multigather_revive_procs(comm, (*sp).ring, local_revive_num, local_revive_procs,
		&(*sp).Revive_procs.num, &(*sp).Revive_procs.procs);
	//if ((*sp).Revive_procs.num > 0) {
	//	printf("total revive num is %d\n", (*sp).Revive_procs.num);
	//}

	remove_revive_procs(sp);
	init_ring(comm_size, my_rank, detector_stage, (*sp).Lagging_procs.num, (*sp).Lagging_procs.procs, &(*sp).ring);
}

//=============================================================================================================================
//(3)激活捡回的进程
//简单地利用链表，捡回后，将进程从链表中删除
//多通信域再考虑其他的实现方式
//=============================================================================================================================

int activate_revive_procs(
	MPI_Comm comm,
	Detector *sp
) {
	if ((*sp).local_re_phead->next == NULL) {
		return 0;
	}

	Revive_LL *pre, *p;
	
	//stage已经发送，仅需要发送lagging数组即可,发送完就删除
	while ((*sp).local_re_phead->next != NULL) {
		p = (*sp).local_re_phead->next;
		pre = (*sp).local_re_phead;

		while (p != NULL) {
			send_procs_template(comm, &p->revive_proc,
				RING_REVIVE_DATA_LAGGING_NUM, RING_REVIVE_DATA_LAGGING_PROCS,
				(*sp).Lagging_procs.num, (*sp).Lagging_procs.procs);
			if (p->revive_proc.comm_stage == FINISH) {
				printf("activate finish\n");
				//如果完成，删除进程
				pre->next = p->next;
				free(p);
				p = pre->next;
			}
			else {
				//如果未完成，指针前移
				pre = p;
				p = p->next;
			}
		}
	}

	return 0;
}





