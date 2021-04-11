#ifndef HEADER_FILE
#define HEADER_FILE
#include "ring_head.h"
#endif

int my_parent(int count, int index) {
	if ((index == 0) || (count == 0))
		return EMPTY;
	else
	{
		if (index % 2 == 0) {
			return (index - 2) / 2;
		}
		else {
			return (index - 1) / 2;
		}
	}
}

int my_lchild(int count, int index) {
	if ((2 * index + 1 > count - 1) || (count == 0)) {
		return EMPTY;
	}
	else {
		return 2 * index + 1;
	}
}

int my_rchild(int count, int index) {
	if ((2 * index + 2 > count - 1) || (count == 0)) {
		return EMPTY;
	}
	else {
		return 2 * index + 2;
	}
}

int tol_tree_bcast(
	int lagging_num,
	int lagging_procs[],
	void *data_recv,
	int count,
	MPI_Datatype datatype,
	int root_proc,
	MPI_Comm comm
) {
	int my_rank, comm_size, i;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &comm_size);

	int p_ind, l_ind, r_ind, my_ind;
	int p_rank, l_rank, r_rank;
	MPI_Request l_req, r_req;
	int l_flag, r_flag;
	int normal_num = comm_size - lagging_num;
	
	if (my_rank == root_proc) {
		my_ind = 0;
	}
	else {
		my_ind = get_comm_index(my_rank, lagging_num, lagging_procs, comm_size);
		if (my_ind == 0) {
			my_ind = get_comm_index(root_proc, lagging_num, lagging_procs, comm_size);
		}
	}
	
	p_ind = my_parent(normal_num, my_ind);
	l_ind = my_lchild(normal_num, my_ind);
	r_ind = my_rchild(normal_num, my_ind);

	if (p_ind != EMPTY) {
		p_rank = get_real_rank(p_ind, lagging_num, lagging_procs, comm_size);
		if (p_ind == 0 && p_rank != root_proc) {
			p_rank = root_proc;
		}
		else if (p_rank == root_proc) {
			p_rank = get_real_rank(0, lagging_num, lagging_procs, comm_size);
		}
	}
	else {
		p_rank = EMPTY;
	}

	if (l_ind != EMPTY) {
		l_rank = get_real_rank(l_ind, lagging_num, lagging_procs, comm_size);
		if (l_rank == root_proc) {
			l_rank = get_real_rank(0, lagging_num, lagging_procs, comm_size);
		}
	}
	else {
		l_rank = EMPTY;
	}

	if (r_ind != EMPTY) {
		r_rank = get_real_rank(r_ind, lagging_num, lagging_procs, comm_size);
		if (r_rank == root_proc) {
			r_rank = get_real_rank(0, lagging_num, lagging_procs, comm_size);
		}
	}
	else {
		r_rank = EMPTY;
	}

	//printf("rank %d, my_ind %d:\n p_ind %d,p_rank %d | l_ind %d,l_rank %d | r_ind %d,r_rank %d\n",
	//	my_rank,my_ind,p_ind,p_rank,l_ind,l_rank,r_ind,r_rank);

	//Ω” ’
	if (p_rank != EMPTY) {
		MPI_Recv(data_recv, count, datatype, p_rank, TREE_BCAST, comm, MPI_STATUSES_IGNORE);
	}
	//∑¢ÀÕ
	if (l_rank != EMPTY) {
		MPI_Issend(data_recv, count, datatype, l_rank, TREE_BCAST, comm, &l_req);
		l_flag = 0;
	}
	else { l_flag = 1; }

	if (r_rank != EMPTY) {
		MPI_Issend(data_recv, count, datatype, r_rank, TREE_BCAST, comm, &r_req);
		r_flag = 0;
	}
	else { r_flag = 1; }

	while ((l_flag && r_flag) != 1) {
		if (l_flag == 0) {
			MPI_Test(&l_req, &l_flag, MPI_STATUSES_IGNORE);
		}

		if (r_flag == 0) {
			MPI_Test(&r_req, &r_flag, MPI_STATUSES_IGNORE);
		}
	}

	return MPI_SUCCESS;
}

