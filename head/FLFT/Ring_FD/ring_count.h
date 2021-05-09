#ifndef HEADER_FILE
#define HEADER_FILE
#include "ring_head.h"
#endif

typedef struct ring_counter {
	int demand;

	int stage;
	int count;
	int sum;
	Ring r;
}Ring_Counter;

void init_ring_counter(Ring_Counter *c) {
	(*c).demand = 0;
	(*c).stage = READY;
	(*c).r.left_proc.rank = -1;
	(*c).r.right_proc.rank = -1;
	(*c).sum = 0;
	(*c).count = 0;
}

void ring_arrived_procs_count(
	MPI_Comm comm,
	int my_rank,
	Ring_Counter *cou
) {
	switch ((*cou).stage)
	{
	case STAGE1: {
		if ((*cou).r.left_proc.rank > my_rank) {
			(*cou).count = 1;
			MPI_Isend(&(*cou).count, 1, MPI_INT, (*cou).r.right_proc.rank, RING_COUNT_REDUCE, comm, &(*cou).r.right_proc.req);
			(*cou).stage = STAGE2;
		}
		else {
			MPI_Iprobe((*cou).r.left_proc.rank, RING_COUNT_REDUCE, comm, &(*cou).r.left_proc.flag, MPI_STATUSES_IGNORE);
			if ((*cou).r.left_proc.flag == 1) {
				MPI_Recv(&(*cou).count, 1, MPI_INT, (*cou).r.left_proc.rank, RING_COUNT_REDUCE, comm, MPI_STATUSES_IGNORE);
				(*cou).count += 1;
				MPI_Isend(&(*cou).count, 1, MPI_INT, (*cou).r.right_proc.rank, RING_COUNT_REDUCE, comm, &(*cou).r.right_proc.req);
				(*cou).stage = STAGE2;
			}
		}
		break;
	}
	case STAGE2: {
		MPI_Test(&(*cou).r.right_proc.req, &(*cou).r.right_proc.flag, MPI_STATUSES_IGNORE);
		if ((*cou).r.right_proc.flag == 1) {
			(*cou).stage = STAGE3;
		}
		break;
	}
	case STAGE3: {
		if ((*cou).r.left_proc.rank > my_rank) {
			MPI_Iprobe((*cou).r.left_proc.rank, RING_COUNT_REDUCE, comm, &(*cou).r.left_proc.flag, MPI_STATUSES_IGNORE);
			if ((*cou).r.left_proc.flag == 1) {
				MPI_Recv(&(*cou).count, 1, MPI_INT, (*cou).r.left_proc.rank, RING_COUNT_REDUCE, comm, MPI_STATUSES_IGNORE);
				(*cou).sum = (*cou).count;

				//printf("rank %d sum is %d\n", my_rank, (*cou).sum);
				MPI_Isend(&(*cou).sum, 1, MPI_INT, (*cou).r.right_proc.rank, RING_COUNT_BCAST, comm, &(*cou).r.right_proc.req);
				(*cou).stage = STAGE4;
			}
		}
		else {
			MPI_Iprobe((*cou).r.left_proc.rank, RING_COUNT_BCAST, comm, &(*cou).r.left_proc.flag, MPI_STATUSES_IGNORE);
			if ((*cou).r.left_proc.flag == 1) {
				MPI_Recv(&(*cou).sum, 1, MPI_INT, (*cou).r.left_proc.rank, RING_COUNT_BCAST, comm, MPI_STATUSES_IGNORE);
				
				//printf("rank %d sum is %d\n", my_rank, (*cou).sum);
				if (my_rank < (*cou).r.right_proc.rank) {
					MPI_Isend(&(*cou).sum, 1, MPI_INT, (*cou).r.right_proc.rank, RING_COUNT_BCAST, comm, &(*cou).r.right_proc.req);
					(*cou).stage = STAGE4;
				}
				else {
					(*cou).stage = FINISH;
				}
			}
		}
		break;
	}
	case STAGE4: {
		MPI_Test(&(*cou).r.right_proc.req, &(*cou).r.right_proc.flag, MPI_STATUSES_IGNORE);
		if ((*cou).r.right_proc.flag == 1) {
			(*cou).stage = FINISH;
		}
		break;
	}
	default:
		break;
	}
}