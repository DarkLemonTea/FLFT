#ifndef HEADER_FILE
#define HEADER_FILE
#include"../fl_head.h"
#endif
#include"../Ring_FD/ring_retrieve.h"

//ͬ���ۼ�revive���̼���
int tree_multigather_revive_procs(
	MPI_Comm comm,
	Ring_Tree rt,
	Last_con_sta lcs,
	P_set ring_revive,
	//output
	P_set *revive
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
		ȱʧ�������ҽ�����Ҫ
	*/
	switch (rt.ide)
	{
	case TREE_ROOT: {
		//���ۼ�����һ�����������Һ��ӵ�revive������Ϣ
		while (1)
		{
			//���ӽ��̵�����,��ȡ����Ϊ����̣��������Ϊ�ҽ���
			data_transmit(comm, rt.tree.left_child, lcs.l_lcs, &lc_tp,
				MISSING_GATHER_DATA_REQUEST, TREE_GATHER_REVIVE_NUM, TREE_GATHER_REVIVE_PROCS,
				&lrecv_set);

			//�Һ��ӽ��̵�����,��ȡ����Ϊ�ҽ��̣��������Ϊ�����
			data_transmit(comm, rt.tree.right_child, lcs.r_lcs, &rc_tp,
				MISSING_GATHER_DATA_REQUEST, TREE_GATHER_REVIVE_NUM, TREE_GATHER_REVIVE_PROCS,
				&rrecv_set);

			if (lc_tp.stage == FINISH && rc_tp.stage == FINISH) {
				lchild_proc.comm_stage = READY; rchild_proc.comm_stage = READY;
				break;
			}
		}

		//���ۼ����ڶ������ϲ�,�õ��ܵ�L����
		merge_tri_procs(ring_revive, lrecv_set, rrecv_set, revive);

		//��������������������������������������������������������������������������������������������������������������������
		//���ಥ����һ��,���͸����Һ��ӽ���		
		while (1)
		{
			if (rt.tree.left_child.rank == EMPTY || lcs.l_lcs == 0) {
				lchild_proc.comm_stage = FINISH;
			}
			else {
				send_procs_template(comm, &lchild_proc, 
					TREE_MULTICAST_REVIVE_NUM, TREE_MULTICAST_REVIVE_PROCS, *revive);
			}

			if (rt.tree.right_child.rank == EMPTY || lcs.r_lcs == 0) {
				rchild_proc.comm_stage = FINISH;
			}
			else {
				send_procs_template(comm, &rchild_proc, 
					TREE_MULTICAST_REVIVE_NUM, TREE_MULTICAST_REVIVE_PROCS, *revive);
			}

			if (lchild_proc.comm_stage == FINISH && rchild_proc.comm_stage == FINISH) {
				break;
			}
		}

		break;
	}
	case TREE_BRANCH: {

		//���ۼ�����һ�����������Һ��ӵ�revive������Ϣ
		while (1)
		{
			//���ӽ��̵�����,��ȡ����Ϊ����̣��������Ϊ�ҽ���
			data_transmit(comm, rt.tree.left_child, lcs.l_lcs, &lc_tp,
				MISSING_GATHER_DATA_REQUEST, TREE_GATHER_REVIVE_NUM, TREE_GATHER_REVIVE_PROCS,
				&lrecv_set);

			//�Һ��ӽ��̵�����,��ȡ����Ϊ�ҽ��̣��������Ϊ�����
			data_transmit(comm, rt.tree.right_child, lcs.r_lcs, &rc_tp,
				MISSING_GATHER_DATA_REQUEST, TREE_GATHER_REVIVE_NUM, TREE_GATHER_REVIVE_PROCS,
				&rrecv_set);

			if (lc_tp.stage == FINISH && rc_tp.stage == FINISH) {
				lchild_proc.comm_stage = READY; rchild_proc.comm_stage = READY;
				break;
			}
		}

		//���ۼ����ڶ������ϲ������ϴ���
		merge_tri_procs(ring_revive, lrecv_set, rrecv_set, &send_set);
		
		while (1)
		{
			if (rt.tree.parent.rank == EMPTY || lcs.p_lcs == 0) {
				break;
			}
			else {
				send_procs_template(comm, &parent_proc, 
					TREE_GATHER_REVIVE_NUM, TREE_GATHER_REVIVE_PROCS, send_set);
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
				MISSING_MULTICAST_DATA_REQUEST, TREE_MULTICAST_REVIVE_NUM, 
				TREE_MULTICAST_REVIVE_PROCS, revive);
			if (p_tp.stage == FINISH) { 
				lchild_proc.comm_stage = READY; rchild_proc.comm_stage = READY;
				break; 
			}
		}

		//���ಥ���ڶ��������͸����Һ���		
		while (1)
		{
			if (rt.tree.left_child.rank == EMPTY || lcs.l_lcs == 0) {
				lchild_proc.comm_stage = FINISH;
			}
			else {
				send_procs_template(comm, &lchild_proc, 
					TREE_MULTICAST_REVIVE_NUM, TREE_MULTICAST_REVIVE_PROCS, *revive);
			}

			if (rt.tree.right_child.rank == EMPTY || lcs.r_lcs == 0) {
				rchild_proc.comm_stage = FINISH;
			}
			else {
				send_procs_template(comm, &rchild_proc, 
					TREE_MULTICAST_REVIVE_NUM, TREE_MULTICAST_REVIVE_PROCS, *revive);
			}

			if (lchild_proc.comm_stage == FINISH && rchild_proc.comm_stage == FINISH) {
				break;
			}
		}

		break;
	}
	case TREE_LEAF: {
		//���ۼ�����һ�������ϴ��ݻ���revive��������
		while (1)
		{
			if (rt.tree.parent.rank == EMPTY || lcs.p_lcs == 0) {
				break;
			}
			else {
				send_procs_template(comm, &parent_proc, 
					TREE_GATHER_REVIVE_NUM, TREE_GATHER_REVIVE_PROCS, ring_revive);
				if (parent_proc.comm_stage == FINISH) {
					parent_proc.comm_stage = READY;
					break;
				}
			}
		}

		//��������������������������������������������������������������������������������������������������������������������
		//���ಥ����һ��,��������˫�׽��̵�L���ϣ�δ���ӳɹ����ҽ�����ȡ
		while (1)
		{
			//����˫�׽������ݣ���ȡ����Ϊ����̣��������Ϊ�ҽ���
			data_transmit(comm, rt.tree.parent, lcs.p_lcs, &p_tp,
				MISSING_MULTICAST_DATA_REQUEST, TREE_MULTICAST_REVIVE_NUM, TREE_MULTICAST_REVIVE_PROCS,
				revive);
			if (p_tp.stage == FINISH) { break; }
		}
		break;
	}
	default:
		break;
	}

	return 0;
}

void rt_remove_revive_procs(
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
			if (in((*sp).Revive_procs, (*sp).Lagging_procs.procs[i])) {
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

//�ۼ���ؽ���
void rt_retrieve_procs(
	int detector_stage,
	MPI_Comm comm,
	double T_retrieve,
	Detector *sp
) {
	int my_rank, comm_size;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &comm_size);
	
	(*sp).Ring_Revive_procs.num = 0;
	(*sp).Ring_Revive_procs.procs = NULL;

	(*sp).Revive_procs.num = 0;
	(*sp).Revive_procs.procs = NULL;

	//��һ������ؽ��̣�ͨ��ͨ��ȷ��
	struct timeval start, end;
	double cost_time = 0;
	gettimeofday(&start, NULL);
	Revive_LL *p, *pre;
	while (cost_time < T_retrieve)
	{
		probe_revive_procs(detector_stage, comm, sp);

		p = (*sp).local_re_phead->next;
		while (p != NULL) {
			rescue_lagging_procs(detector_stage, comm, &(p->revive_proc));
			p = p->next;
		}
		gettimeofday(&end, NULL);
		cost_time = 1000.0 * (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000.0;
	}

	//�ڶ�����ɸѡ����صı��ؽ��̣����ʧ�ܵ�ɾ��
	P_set local_revive;
	init_p_set(&local_revive);

	//��һ�α���������
	pre = (*sp).local_re_phead;
	p = (*sp).local_re_phead->next;
	while (p != NULL) {
		if (p->revive_proc.comm_stage == FINISH) {
			local_revive.num += 1;
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
	local_revive.procs = (int*)calloc(local_revive.num, sizeof(int));

	//�ڶ��α���������
	pre = (*sp).local_re_phead;
	p = (*sp).local_re_phead->next;
	int i = 0;
	while (p != NULL) {
		local_revive.procs[i] = p->revive_proc.rank;
		i += 1;
		pre = p;
		p = p->next;
	}

	//printf("rank %d local revive num is %d\n", my_rank, local_revive.num);

	ring_multigather_revive_procs(comm, (*sp).rt.ring, local_revive, &(*sp).Ring_Revive_procs);

	//printf("rank %d ring revive num is %d\n", my_rank, (*sp).Ring_Revive_procs.num);
	
	tree_multigather_revive_procs(comm, (*sp).rt, (*sp).rt.lcs, 
		(*sp).Ring_Revive_procs, &(*sp).Revive_procs);

	/*printf("rank %d rln is %d,ln is %d before procs removed\n",
		my_rank, (*sp).Ring_lagging_procs.num, (*sp).Lagging_procs.num);*/
	
	/*printf("rank %d rrn is %d,rn is %d before procs removed\n",
		my_rank, (*sp).Ring_Revive_procs.num, (*sp).Revive_procs.num);*/

	remove_sub_ring_revive_procs(sp);

	rt_remove_revive_procs(sp);

	/*printf("rank %d rln is %d,ln is %d after procs removed\n",
		my_rank, (*sp).Ring_lagging_procs.num, (*sp).Lagging_procs.num);*/

	/*printf("rank %d rrn is %d,rn is %d after procs removed\n",
		my_rank, (*sp).Ring_Revive_procs.num, (*sp).Revive_procs.num);*/

	init_ring_tree(comm_size, my_rank, (*sp).rt.ring_scale, detector_stage,
		(*sp).Ring_lagging_procs, &(*sp).rt);

	/*printf("rank %d ring left is %d,ring right is %d\n",
		my_rank, (*sp).rt.ring.left_proc.rank, (*sp).rt.ring.right_proc.rank);*/
}






