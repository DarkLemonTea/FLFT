#include <string.h>

#ifndef HEADER_FILE
#define HEADER_FILE
#include "../Ring_FD/ring_FD.h"
#endif

//int bruck_allgather(
//	void *sendbuf,
//	int sendcount,
//	MPI_Datatype sendtype,
//	void *recvbuf,
//	int recvcount,
//	MPI_Datatype recvtype,
//	MPI_Comm comm
//){
//	int my_rank, comm_size, sendto, recvfrom, distance, blockcount;
//	MPI_Aint slb, rlb, sext, rext;
//	char *tmpsend = NULL, *tmprecv = NULL;
//
//	int send_type_size, recv_type_size;
//	MPI_Type_size(sendtype, &send_type_size);
//	MPI_Type_size(recvtype, &recv_type_size);
//
//	MPI_Comm_rank(comm, &my_rank);
//	MPI_Comm_size(comm, &comm_size);
//	
//	MPI_Type_get_extent(sendtype, &slb, &sext);
//	MPI_Type_get_extent(recvtype, &rlb, &rext);
//
//	tmprecv = (char*)recvbuf;
//	if (MPI_IN_PLACE != sendbuf) {
//		tmpsend = (char*)sendbuf;
//		memcpy(tmprecv, tmpsend, recvcount * recv_type_size);
//	}
//	else if (0 != my_rank) { 
//		tmpsend = ((char*)recvbuf) + (MPI_Aint)my_rank * (MPI_Aint)recvcount * rext;
//		memcpy(tmprecv, tmpsend, recvcount * recv_type_size);
//	}
//	
//	blockcount = 1;
//	tmpsend = (char*)recvbuf;
//	for (distance = 1; distance < comm_size; distance <<= 1) {
//
//		recvfrom = (my_rank + distance) % comm_size;
//		sendto = (my_rank - distance + comm_size) % comm_size;
//
//		tmprecv = tmpsend + (MPI_Aint)distance * (MPI_Aint)recvcount * rext;
//
//		if (distance <= (comm_size >> 1)) {
//			blockcount = distance;
//		}
//		else {
//			blockcount = comm_size - distance;
//		}
//
//		/* Sendreceive */
//		MPI_Sendrecv(tmpsend, blockcount * recvcount, recvtype, sendto, BRUCK_ALLGATHER,
//			tmprecv, blockcount * recvcount, recvtype, recvfrom, BRUCK_ALLGATHER,
//			comm, MPI_STATUSES_IGNORE);
//	}
//
//	if (0 != my_rank) {
//		char *free_buf = NULL, *shift_buf = NULL;
//		MPI_Aint extent, true_extent, lb, true_lb;
//		MPI_Aint span, gap = 0;
//
//		MPI_Type_get_extent(recvtype, &lb, &extent);
//		MPI_Type_get_extent(recvtype, &true_lb, &true_extent);
//		
//		gap = true_lb;
//		span = true_extent + (recvcount - 1) * extent;
//
//		free_buf = (char*)calloc(span, sizeof(char));
//		shift_buf = free_buf - gap;
//
//		/* 1. copy blocks [0 .. (size - rank - 1)] from rbuf to shift recver */
//		memcpy(shift_buf, recvbuf, 
//			recv_type_size * ((MPI_Aint)(comm_size - my_rank) * (MPI_Aint)recvcount));
//
//		/* 2. move blocks [(size - rank) .. size] from rbuf to the begining of rbuf */
//		tmpsend = (char*)recvbuf + (MPI_Aint)(comm_size - my_rank) * (MPI_Aint)recvcount * rext;
//		memcpy(recvbuf, tmpsend,
//			recv_type_size * ((MPI_Aint)my_rank * (MPI_Aint)recvcount));
//
//		/* 3. copy blocks from shift recver back to rbuf starting at block [rank]. */
//		tmprecv = (char*)recvbuf + (MPI_Aint)my_rank * (MPI_Aint)recvcount * rext;
//		memcpy(tmprecv, shift_buf,
//			recv_type_size * ((MPI_Aint)(comm_size - my_rank) * (MPI_Aint)recvcount));
//
//		free(free_buf);
//	}
//
//	return 0;
//}

int get_handle_num(
	P_set Lagging_procs,
	int my_rank,
	int comm_size
) {
	if ( Lagging_procs.num == 0) return 0;

	int ind = comm_size, i, handle_num = 0;
	for (i = 0; i < Lagging_procs.num; i++) {
		if (Lagging_procs.procs[i] == (my_rank - 1 + comm_size) % comm_size) {
			ind = i;
			handle_num += 1;
			break;
		}
	}
	if (ind < comm_size) {
		while (Lagging_procs.procs[ind] - 1 == Lagging_procs.procs[ind - 1])
		{
			ind -= 1;
			handle_num += 1;
		}
		if (ind == 0 && Lagging_procs.procs[0] == 0 && Lagging_procs.procs[Lagging_procs.num - 1] == comm_size - 1) {
			handle_num += 1;
			ind = Lagging_procs.num - 1;
			while (Lagging_procs.procs[ind] - 1 == Lagging_procs.procs[ind - 1])
			{
				ind -= 1;
				handle_num += 1;
			}
		}
	}
	
	return handle_num;
}

int tol_bruck_multigather(
	P_set Lagging_procs,
	int diff, //偏移数，即每个进程负责故障进程数目
	void *sendbuf, 
	MPI_Datatype sendtype,
	void *recvbuf,
	MPI_Datatype recvtype,
	int datacount, //每组数据数目
	MPI_Comm comm
) {
	//排序
	QuickSort(Lagging_procs.procs, 0, Lagging_procs.num - 1);

	int my_rank, comm_size, distance;
	MPI_Request req;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &comm_size);

	int send_ind, send_to, send_num;
	int recv_ind, recv_from, recv_num;
	MPI_Aint slb, rlb, sext, rext;
	int n = comm_size - Lagging_procs.num;	//参与多播的数量n

	char *tmpsend = NULL, *tmprecv = NULL;
	int my_ind = get_comm_index(my_rank, Lagging_procs.num, Lagging_procs.procs, comm_size);

	int send_type_size, recv_type_size;
	MPI_Type_size(sendtype, &send_type_size);
	MPI_Type_size(recvtype, &recv_type_size);

	MPI_Type_get_extent(sendtype, &slb, &sext);
	MPI_Type_get_extent(recvtype, &rlb, &rext);

	//不支持MPI_IN_PLACE参数
	tmprecv = (char*)recvbuf;
	tmpsend = (char*)sendbuf;
	send_num = diff + 1;

	memcpy(tmprecv, tmpsend, send_num * datacount * recv_type_size);

	tmpsend = (char*)recvbuf;
	for (distance = 1; distance < n; distance <<= 1) {

		recv_ind = (my_ind + distance) % n;
		send_ind = (my_ind + n - distance) % n;
		recv_from = get_real_rank(recv_ind, Lagging_procs.num, Lagging_procs.procs, comm_size);
		send_to = get_real_rank(send_ind, Lagging_procs.num, Lagging_procs.procs, comm_size);

		/*printf("rank %d, my_ind %d: dis %d,recv_ind is %d,recv_from is %d\n", 
			my_rank, my_ind, distance, recv_ind, recv_from);*/
		
		/*printf("rank %d, my_ind %d: dis %d,send_ind is %d,send_to is %d\n", 
			my_rank, my_ind, distance, send_ind, send_to);*/	

		tmprecv = tmpsend + (MPI_Aint)send_num * (MPI_Aint)datacount * rext;
		if (distance * 2 <= n) {
			MPI_Irecv(&recv_num, 1, MPI_INT, recv_from, BRUCK_RECVNUM, comm, &req);
			MPI_Send(&send_num, 1, MPI_INT, send_to, BRUCK_RECVNUM, comm);
			MPI_Wait(&req, MPI_STATUSES_IGNORE);		
			//printf("rank %d: dis %d, send_num %d, recv_num is %d\n", my_rank, distance, send_num, recv_num);
			
			MPI_Irecv(tmprecv, recv_num * datacount, recvtype, recv_from, BRUCK_MULTIGATHER, comm, &req);
			MPI_Send(tmpsend, send_num * datacount, recvtype, send_to, BRUCK_MULTIGATHER, comm);
			MPI_Wait(&req, MPI_STATUSES_IGNORE);

			send_num += recv_num;		
		}
		else {
			//printf("rank %d send num is %d\n", my_rank, send_num);
			recv_num = comm_size - send_num;
			MPI_Irecv(&send_num, 1, MPI_INT, send_to, BRUCK_RECVNUM, comm, &req);
			MPI_Send(&recv_num, 1, MPI_INT, recv_from, BRUCK_RECVNUM, comm);
			MPI_Wait(&req, MPI_STATUSES_IGNORE);
			//printf("rank %d: dis %d, send_num %d, recv_num is %d\n", my_rank, distance, send_num, recv_num);

			MPI_Irecv(tmprecv, recv_num * datacount, recvtype, recv_from, BRUCK_MULTIGATHER, comm, &req);
			MPI_Send(tmpsend, send_num * datacount, recvtype, send_to, BRUCK_MULTIGATHER, comm);
			MPI_Wait(&req, MPI_STATUSES_IGNORE);
		}	
	}

	//调整顺序
	int first_count = (comm_size - my_rank + diff) % comm_size;
	char *free_buf = NULL, *shift_buf = NULL;
	MPI_Aint extent, true_extent, lb, true_lb;
	MPI_Aint span, gap = 0;

	MPI_Type_get_extent(recvtype, &lb, &extent);
	MPI_Type_get_extent(recvtype, &true_lb, &true_extent);

	gap = true_lb;
	span = (true_extent + (first_count * datacount - 1) * extent);

	free_buf = (char*)calloc(span, sizeof(char));
	shift_buf = free_buf - gap;

	/* 1. copy blocks [0 .. (size - rank - 1)] from rbuf to shift recver */
	memcpy(shift_buf, recvbuf,
		recv_type_size * ((MPI_Aint)(first_count) * (MPI_Aint)datacount));

	/* 2. move blocks [(size - rank) .. size] from rbuf to the begining of rbuf */
	tmpsend = (char*)recvbuf + (MPI_Aint)(first_count) * (MPI_Aint)datacount * (MPI_Aint)recv_type_size;
	memcpy(recvbuf, tmpsend,
		recv_type_size * ((MPI_Aint)(comm_size - first_count) * (MPI_Aint)datacount));

	/* 3. copy blocks from shift recver back to rbuf starting at block [rank]. */
	tmprecv = (char*)recvbuf + (MPI_Aint)(comm_size - first_count) * (MPI_Aint)datacount * (MPI_Aint)recv_type_size;
	memcpy(tmprecv, shift_buf,
		recv_type_size * ((MPI_Aint)(first_count) * (MPI_Aint)datacount));

	free(free_buf);

	return 0;
}
