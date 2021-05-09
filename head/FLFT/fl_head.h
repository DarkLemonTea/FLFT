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

#define TYPE_NONE 9			//���ݴ�
#define TYPE_RING 10		//������
#define TYPE_RING_TREE 11		//��-������
#define TYPE_RING_BUTTERFLY 12		//��-��������

//�ȴ�ʱ������б���
typedef struct FD_variable {
	double release_rate;
	int max_fail_num;
	double T_max;
	double T_wait;
	double T_retrieve;
}FD_var;

void init_fd_var(
	double release_rate, //������
	int max_fail_num, //�����Ͻ�������
	double T_max,        //wall time ������ʱ�䣬���ش���
	double T_wait,       //wait time ������ʱ�䣬���Ҵﵽ�����ʣ���ɹ��ϼ��
	double T_retrieve,   //���̼�صĵȴ�ʱ��
	FD_var *fd
) {
	//���뵥λ���룬ת��Ϊ����
	(*fd).T_max = T_max * 1000;
	(*fd).T_wait = T_wait * 1000;
	(*fd).T_retrieve = T_retrieve * 1000;
	(*fd).release_rate = release_rate;
	(*fd).max_fail_num = max_fail_num;
}

//=============================================================================================================================
//ͨ�õ�detector������ͨ�Žṹ��ͬѡ���ʼ����ʽ
//=============================================================================================================================

typedef struct Set_pointers {
	/*ͨ�ýṹ*/
	int type;					 //��������
	int current_stage;			 //��ǰ�׶�
	P_set Lagging_procs;         //�ٵ������б�
	Revive_LL *local_re_phead;   //�洢�ָ����̵���ʱ����	
	P_set Revive_procs;          //�ָ������б�

	P_set Ring_lagging_procs;    //ͬ���ڳٵ������б����ڿ��ٳ�ʼ��
	P_set Ring_Revive_procs;     //ͬ���ڻָ������б�

	//��
	Ring ring;
	
	//��-��
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

	//��ʼ��
	(*sp).current_stage = d_stage;
	(*sp).type = type;

	init_p_set(&(*sp).Lagging_procs);
	init_p_set(&(*sp).Revive_procs);

	//����һ��ͷ��㣬�������κ�����
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
//���revive���̵�����
//=============================================================================================================================

//��ͷ��㣨ͷ��㲻�洢���ݣ���ͷ�巨
void creat_and_init(Revive_LL *head, int rank, int stage) {
	Revive_LL *node;
	node = (Revive_LL*)malloc(sizeof(Revive_LL));
	node->revive_proc.rank = rank;
	node->revive_proc.send = stage;
	node->revive_proc.comm_stage = STAGE1;

	node->next = head->next;
	head->next = node;
}

//����
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
//������Ϣ
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
//�����Ӳ���
//===========================================================================================================================
//���ѡһ��ͬ������
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

//ͳ�Ƴٵ����̣�����������Ĳ��֣�
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

//��lagging�����л�ȡͬ���ڽ���
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
		//ͬ��ͬ��
		if (my_rank / ring_scale == p.procs[i] / ring_scale) {
			(*rp).num += 1;
		}
	}
	(*rp).procs = (int*)calloc((*rp).num, sizeof(int));
	for (i = 0; i < p.num; i++) {
		//ͬ��ͬ��
		if (my_rank / ring_scale == p.procs[i] / ring_scale) {
			(*rp).procs[ind] = p.procs[i];
			ind += 1;
		}
	}
}

//===========================================================================================================================
//�������������ݴ��䡢����
//===========================================================================================================================
typedef struct transmit_procs {
	int stage;            //����׶�
	Comm_proc dir_proc;  //����������Ķ��󣬳�ʼ����ͨ�ű���
	Comm_proc req_proc;   //��ȡ����ȱʧʱ�Ҹö�����Ҫ
	Comm_proc give_proc;  //������󣬽���ȡ�����ݺ󣬿��Խ���ȱʧ����
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

//���ݴ���
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
			(*tp).req_proc.send = 0; //����Ҫ�ȱʧ����
			(*tp).stage = RECV_STAGE;
		}
		else if (con_status == 0) {
			(*tp).req_proc.send = 1; //��ʾ���ȱʧ��������
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
		//����ȡ�ŷ����ź�
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

		//����ո��跽�ź�
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

//��������
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
			(*tp).req_proc.send = 0; //����Ҫ�ȱʧ����
			(*tp).stage = RECV_STAGE;
		}
		else if (con_status == 0) {
			(*tp).req_proc.send = 1; //��ʾ���ȱʧ��������
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
		//����ȡ�ŷ����ź�
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

		//����ո��跽�ź�
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