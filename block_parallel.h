#include<stdlib.h>
#include<string.h>
#include<stddef.h>
#include<mpi.h>
#include "AnyArray.h"

#define BLOCK_PARALLEL_REDUCE 50000
#define BLOCK_PARALLEL_BCAST 50001

void free_reqs(MPI_Request *reqs, int num) {
	if (reqs == NULL) {
		return;
	}

	int i = 0;
	for (i = 0; i < num; i++){
		if (reqs[i] != MPI_REQUEST_NULL) {
			MPI_Request_free(&reqs[i]);
		}
	}
}

MPI_Aint datatype_span(MPI_Datatype dt,int count,MPI_Aint *gap) {
	int size;
	MPI_Aint lb, t_lb, extent, t_extent;
	MPI_Type_size(dt, &size);
	MPI_Type_get_extent(dt, &lb, &extent);
	MPI_Type_get_true_extent(dt, &t_lb, &t_extent);
	if (count == 0 || size == 0) {
		*gap = 0;
		return 0;
	}
	*gap = t_lb;
	return t_extent + (count - 1)*extent;
}

typedef struct tree_mem{
	int rank;
}Tree_mem;

typedef struct coll_tree {
	int root;
	Tree_mem parent;
	Tree_mem lchild;
	Tree_mem rchild;
	int nextsize;		//孩子个数
}Coll_Tree;

//用于分块集合通信的树结构，n*n=comm_size
void creat_cube_tree(
	int my_rank,
	int n,
	char row_or_col,
	int root,
	Coll_Tree *ct
) {
	(*ct).nextsize = 0;
	(*ct).root = root;

	int row_ind = my_rank / n;
	int col_ind = my_rank % n;

	int i = 0, index;
	int *procs;
	procs = (int*)calloc(n, sizeof(int));

	if (row_or_col == 'r') {
		for (i = 0; i < n; i++) {
			procs[i] = row_ind * n + i;
		}
		//交换root位置
		procs[0] = root; procs[root%n] = row_ind * n;
		if (my_rank == root) {
			index = 0;
		}
		else {
			if (col_ind == 0) { index = root % n; }
			else { index = my_rank % n; }
		}
		//根据坐标，读取对应树的值
		//树方向：
		if (index == 0) { (*ct).parent.rank = -1; }
		else { (*ct).parent.rank = procs[(index - 1) / 2]; }
		//左孩子
		if (index * 2 + 1 >= n) { (*ct).lchild.rank = -1; }
		else { (*ct).lchild.rank = procs[index * 2 + 1]; }
		//右孩子
		if (index * 2 + 2 >= n) { (*ct).rchild.rank = -1; }
		else { (*ct).rchild.rank = procs[index * 2 + 2]; }
		free(procs);
	}
	else if (row_or_col = 'c') {
		for (i = 0; i < n; i++) {
			procs[i] = i * n + col_ind;
		}
		//交换root位置
		procs[0] = root; procs[root/n] = col_ind;
		if (my_rank == root) {
			index = 0;
		}
		else {
			if (row_ind == 0) { index = root / n; }
			else { index = my_rank / n; }
		}
		//根据坐标，读取对应树的值
		//树方向：
		if (index == 0) { (*ct).parent.rank = -1; }
		else { (*ct).parent.rank = procs[(index - 1) / 2]; }
		//左孩子
		if (index * 2 + 1 >= n) { (*ct).lchild.rank = -1; }
		else { (*ct).lchild.rank = procs[index * 2 + 1]; }
		//右孩子
		if (index * 2 + 2 >= n) { (*ct).rchild.rank = -1; }
		else { (*ct).rchild.rank = procs[index * 2 + 2]; }
		free(procs);
	}

	if ((*ct).lchild.rank > -1) {
		(*ct).nextsize += 1;
	}
	if ((*ct).rchild.rank > -1) {
		(*ct).nextsize += 1;
	}
}

int Reduce_2D(
	void* sendbuf,
	void* recvbuf,
	int original_count,
	MPI_Datatype datatype,
	MPI_Op op,
	Coll_Tree ct,
	MPI_Comm comm,
	MPI_Aint count_by_segment,
	int max_outstanding_reqs
) {
	int my_rank, comm_size;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &comm_size);

	//printf("row rank %d root is %d;r_ct p is %d, l is %d, r is %d\n",
	//	my_rank, ct.root, ct.parent.rank, ct.lchild.rank, ct.rchild.rank);

	char *inbuf[2] = { NULL, NULL }, *inbuf_free[2] = { NULL, NULL };
	char *accumbuf = NULL, *accumbuf_free = NULL;
	char *local_op_buffer = NULL, *sendtmpbuf = NULL;
	MPI_Aint extent, size, gap = 0, segment_increment;
	MPI_Request *sreq = NULL, reqs[2] = { MPI_REQUEST_NULL, MPI_REQUEST_NULL };
	int num_segments, line, ret, segindex, i;
	int recvcount, prevcount, inbi;
	int type_size;

	int child[2] = { ct.lchild.rank, ct.rchild.rank };

	/**
	 * Determine number of segments and number of elements
	 * sent per operation
	 */
	MPI_Type_size(datatype, &type_size);
	MPI_Type_extent(datatype, &extent);
	num_segments = (int)(((size_t)original_count + (size_t)count_by_segment - (size_t)1) / (size_t)count_by_segment);
	segment_increment = (MPI_Aint)count_by_segment * extent;

	sendtmpbuf = (char*)sendbuf;
	if (sendbuf == MPI_IN_PLACE) {
		sendtmpbuf = (char *)recvbuf;
	}

	sendtmpbuf = (char*)sendbuf;

	if (ct.nextsize > 0) {
		ptrdiff_t real_segment_size;

		/* handle non existant recv buffer (i.e. its NULL) and
		   protect the recv buffer on non-root nodes */
		accumbuf = (char*)recvbuf;
		if ((NULL == accumbuf) || (ct.root != my_rank)) {
			/* Allocate temporary accumulator buffer. */
			size = datatype_span(datatype, original_count, &gap);
			accumbuf_free = (char*)malloc(size);
			if (accumbuf_free == NULL) {
				line = __LINE__; ret = -1; goto error_hndl;
			}
			accumbuf = accumbuf_free - gap;
		}

		/* If this is a non-commutative operation we must copy
		   sendbuf to the accumbuf, in order to simplfy the loops */
		   //将sendbuf复制到accumbuf，简化循环
		memcpy(accumbuf, sendtmpbuf, original_count * type_size);

		real_segment_size = datatype_span(datatype, count_by_segment, &gap);
		inbuf_free[0] = (char*)malloc(real_segment_size);
		if (inbuf_free[0] == NULL) {
			line = __LINE__; ret = -1; goto error_hndl;
		}
		// 为传入的段分配两个缓冲区
		inbuf[0] = inbuf_free[0] - gap;
		// 要用到交替通信，分配第二个缓冲区
		if ((num_segments > 1) || (ct.nextsize > 1)) {
			inbuf_free[1] = (char*)malloc(real_segment_size);
			if (inbuf_free[1] == NULL) {
				line = __LINE__; ret = -1; goto error_hndl;
			}
			inbuf[1] = inbuf_free[1] - gap;
		}

		/* reset input buffer index and receive count */
		inbi = 0;
		recvcount = 0;
		/* for each segment */
		for (segindex = 0; segindex <= num_segments; segindex++) {
			prevcount = recvcount;
			/* recvcount - number of elements in current segment */
			recvcount = count_by_segment;
			if (segindex == (num_segments - 1))
				recvcount = original_count - (ptrdiff_t)count_by_segment * (ptrdiff_t)segindex;

			/* for each child */
			for (i = 0; i < ct.nextsize; i++) {
				/**
				 * We try to overlap communication:
				 * either with next segment or with the next child
				 */
				 /* post irecv for current segindex on current child */
				if (segindex < num_segments) {
					void* local_recvbuf = inbuf[inbi];
					if (0 == i) {
						/* for the first step (1st child per segment) and
						 * commutative operations we might be able to irecv
						 * directly into the accumulate buffer so that we can
						 * reduce(op) this with our sendbuf in one step as
						 * ompi_op_reduce only has two buffer pointers,
						 * this avoids an extra memory copy.
						 *
						 * BUT if the operation is non-commutative or
						 * we are root and are USING MPI_IN_PLACE this is wrong!
						 */
						if (my_rank == ct.root) {
							local_recvbuf = 
								accumbuf + (ptrdiff_t)segindex * (ptrdiff_t)segment_increment;
						}
					}

					ret = MPI_Irecv(local_recvbuf, recvcount, datatype, child[i],
						BLOCK_PARALLEL_REDUCE, comm, &reqs[inbi]);
					if (ret != MPI_SUCCESS) { line = __LINE__; goto error_hndl; }
				}
				/* wait for previous req to complete, if any.
				   if there are no requests reqs[inbi ^1] will be
				   MPI_REQUEST_NULL. */
				   /* wait on data from last child for previous segment */
				ret = MPI_Wait(&reqs[inbi ^ 1], MPI_STATUSES_IGNORE);
				if (ret != MPI_SUCCESS) { line = __LINE__; goto error_hndl; }
				
				local_op_buffer = inbuf[inbi ^ 1];
				if (i > 0) {
					/* our first operation is to combine our own [sendbuf] data
					 * with the data we recvd from down stream (but only
					 * the operation is commutative and if we are not root and
					 * not using MPI_IN_PLACE)
					 */
					if (1 == i) {
						if (my_rank == ct.root) {
							local_op_buffer = 
								sendtmpbuf + (ptrdiff_t)segindex * (ptrdiff_t)segment_increment;
						}
					}
					local_reduce(op, local_op_buffer, 
						accumbuf + (ptrdiff_t)segindex * (ptrdiff_t)segment_increment,
						recvcount, datatype);

				}
				else if (segindex > 0) {
					void* accumulator = 
						accumbuf + (ptrdiff_t)(segindex - 1) * (ptrdiff_t)segment_increment;
					if (ct.nextsize <= 1) {
						if (my_rank == ct.root) {
							local_op_buffer = 
								sendtmpbuf + (ptrdiff_t)(segindex - 1) * (ptrdiff_t)segment_increment;
						}
					}
					local_reduce(op, local_op_buffer, accumulator, 
						prevcount, datatype);

					/* all reduced on available data this step (i) complete,
					 * pass to the next process unless you are the root.
					 */
					if (my_rank != ct.root) {
						/* send combined/accumulated data to parent */
						ret = MPI_Send(accumulator, prevcount, datatype, 
							ct.parent.rank, BLOCK_PARALLEL_REDUCE, comm);
						if (ret != MPI_SUCCESS) {
							line = __LINE__; goto error_hndl;
						}
					}

					/* we stop when segindex = number of segments
					   (i.e. we do num_segment+1 steps for pipelining */
					if (segindex == num_segments) break;
				}

				/* update input buffer index */
				inbi = inbi ^ 1;
			} /* end of for each child */
		} /* end of for each segment */

		/* clean up */
		if (inbuf_free[0] != NULL) free(inbuf_free[0]);
		if (inbuf_free[1] != NULL) free(inbuf_free[1]);
		if (accumbuf_free != NULL) free(accumbuf_free);
	}
	else {

	/* If the number of segments is less than a maximum number of oustanding
	   requests or there is no limit on the maximum number of outstanding
	   requests, we send data to the parent using blocking send */
	   if ((0 == max_outstanding_reqs) ||
		   (num_segments <= max_outstanding_reqs)) {

		   segindex = 0;
		   while (original_count > 0) {
			   if (original_count < count_by_segment) {
				   count_by_segment = original_count;
			   }
			   ret = MPI_Send((char*)sendbuf +
				   (ptrdiff_t)segindex * (ptrdiff_t)segment_increment,
				   count_by_segment, datatype,
				   ct.parent.rank, BLOCK_PARALLEL_REDUCE, comm);
			   if (ret != MPI_SUCCESS) { line = __LINE__; goto error_hndl; }
			   segindex++;
			   original_count -= count_by_segment;
		   }
	   }

	   /* Otherwise, introduce flow control:
		  - post max_outstanding_reqs non-blocking synchronous send,
		  - for remaining segments
		  - wait for a ssend to complete, and post the next one.
		  - wait for all outstanding sends to complete.
	   */
	   else {
		   int creq = 0;
		   sreq = (MPI_Request*)calloc(max_outstanding_reqs, sizeof(MPI_Request));
		   if (NULL == sreq) { line = __LINE__; ret = -1; goto error_hndl; }

		   /* post first group of requests */
		   for (segindex = 0; segindex < max_outstanding_reqs; segindex++) {
			   ret = MPI_Isend((char*)sendbuf + (ptrdiff_t)segindex * (ptrdiff_t)segment_increment,
				   count_by_segment, datatype, ct.parent.rank,
				   BLOCK_PARALLEL_REDUCE, comm, &sreq[segindex]);
				   
			   if (ret != MPI_SUCCESS) { line = __LINE__; goto error_hndl; }
			   original_count -= count_by_segment;
		   }

		   creq = 0;
		   while (original_count > 0) {
			   /* wait on a posted request to complete */
			   ret = MPI_Wait(&sreq[creq], MPI_STATUS_IGNORE);
			   if (ret != MPI_SUCCESS) { line = __LINE__; goto error_hndl; }

			   if (original_count < count_by_segment) {
				   count_by_segment = original_count;
			   }
			   ret = MPI_Isend((char*)sendbuf + (ptrdiff_t)segindex * (ptrdiff_t)segment_increment,
				   count_by_segment, datatype, ct.parent.rank,
				   BLOCK_PARALLEL_REDUCE, comm, &sreq[creq]);
			   if (ret != MPI_SUCCESS) { line = __LINE__; goto error_hndl; }
			   creq = (creq + 1) % max_outstanding_reqs;
			   segindex++;
			   original_count -= count_by_segment;
		   }

		   /* Wait on the remaining request to complete */
		   ret = MPI_Waitall(max_outstanding_reqs, sreq, MPI_STATUSES_IGNORE);
		   if (ret != MPI_SUCCESS) { line = __LINE__; goto error_hndl; }
	   }
	}

	return MPI_SUCCESS;

error_hndl:  // 错误处理程序 error handler
	printf("ERROR_HNDL: node %d file %s line %d error %d\n", my_rank, __FILE__, line, ret);
	(void)line;  // silence compiler warning
	if (inbuf[0] != NULL) free(inbuf[0]);
	if (inbuf[1] != NULL) free(inbuf[1]);
	if (accumbuf != NULL) free(accumbuf);
	if (NULL != sreq) {
		free_reqs(sreq, max_outstanding_reqs);
		free(sreq);
	}
	
	return ret;
}

int Bcast_2D(
	void* buffer,
	int original_count,
	MPI_Datatype datatype,
	Coll_Tree ct,
	MPI_Comm comm,
	int count_by_segment
) {
	int my_rank, comm_size;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &comm_size);

	//printf("col rank %d root is %d;b_ct p is %d, l is %d, r is %d\n",
	//	my_rank, ct.root, ct.parent.rank, ct.lchild.rank, ct.rchild.rank);

	int err = 0, line, i, segindex, req_index;
	int num_segments; /* Number of segments */
	int sendcount;    /* number of elements sent in this segment */
	size_t realsegsize;
	MPI_Aint extent;
	int type_size;
	char *tmpbuf;	
	
	MPI_Request *send_reqs = NULL, recv_reqs[2] = { MPI_REQUEST_NULL, MPI_REQUEST_NULL };

	int child[2] = { ct.lchild.rank,ct.rchild.rank };

	MPI_Type_size(datatype, &type_size);
	MPI_Type_extent(datatype, &extent);
	num_segments = (int)(((size_t)original_count + (size_t)count_by_segment - (size_t)1) / (size_t)count_by_segment);
	realsegsize = (ptrdiff_t)count_by_segment * extent;

	tmpbuf = (char*)buffer;

	if (ct.nextsize > 0) {
		send_reqs = (MPI_Request*)malloc(ct.nextsize * sizeof(MPI_Request));
	}

	/* Root code */
	if (my_rank == ct.root) {
		/*
		   For each segment:
		   - send segment to all children.
		   The last segment may have less elements than other segments.
		*/
		sendcount = count_by_segment;
		for (segindex = 0; segindex < num_segments; segindex++) {
			if (segindex == (num_segments - 1)) {
				sendcount = original_count - segindex * count_by_segment;
			}
			for (i = 0; i < ct.nextsize; i++) {
				err = MPI_Isend(tmpbuf, sendcount, datatype, child[i],
					BLOCK_PARALLEL_BCAST, comm, &send_reqs[i]);
					
				if (err != MPI_SUCCESS) { line = __LINE__; goto error_hndl; }
			}

			/* complete the sends before starting the next sends */
			err = MPI_Waitall(ct.nextsize, send_reqs, MPI_STATUSES_IGNORE);
			if (err != MPI_SUCCESS) { line = __LINE__; goto error_hndl; }

			/* update tmp buffer */
			tmpbuf += realsegsize;
		}
	}
	else if (ct.nextsize > 0) {
		/*
		   Create the pipeline.
		   1) Post the first receive
		   2) For segments 1 .. num_segments
		   - post new receive
		   - wait on the previous receive to complete
		   - send this data to children
		   3) Wait on the last segment
		   4) Compute number of elements in last segment.
		   5) Send the last segment to children
		*/
		req_index = 0;
		err = MPI_Irecv(tmpbuf, count_by_segment, datatype,
			ct.parent.rank, BLOCK_PARALLEL_BCAST, comm, &recv_reqs[req_index]);
			
		if (err != MPI_SUCCESS) { line = __LINE__; goto error_hndl; }

		for (segindex = 1; segindex < num_segments; segindex++) {

			req_index = req_index ^ 0x1;

			/* post new irecv */
			err = MPI_Irecv(tmpbuf + realsegsize, count_by_segment,
				datatype, ct.parent.rank, BLOCK_PARALLEL_BCAST,
				comm, &recv_reqs[req_index]);
			if (err != MPI_SUCCESS) { line = __LINE__; goto error_hndl; }

			/* wait for and forward the previous segment to children */
			err = MPI_Wait(&recv_reqs[req_index ^ 0x1],
				MPI_STATUS_IGNORE);
			if (err != MPI_SUCCESS) { line = __LINE__; goto error_hndl; }

			for (i = 0; i < ct.nextsize; i++) {
				err = MPI_Isend(tmpbuf, count_by_segment, datatype,
					child[i], BLOCK_PARALLEL_BCAST, comm, &send_reqs[i]);
				if (err != MPI_SUCCESS) { line = __LINE__; goto error_hndl; }
			}

			/* complete the sends before starting the next iteration */
			err = MPI_Waitall(ct.nextsize, send_reqs,
				MPI_STATUSES_IGNORE);
			if (err != MPI_SUCCESS) { line = __LINE__; goto error_hndl; }

			/* Update the receive buffer */
			tmpbuf += realsegsize;
		}

		/* Process the last segment */
		err = MPI_Wait(&recv_reqs[req_index], MPI_STATUS_IGNORE);
		if (err != MPI_SUCCESS) { line = __LINE__; goto error_hndl; }
		
		sendcount = original_count - (ptrdiff_t)(num_segments - 1) * count_by_segment;
		for (i = 0; i < ct.nextsize; i++) {
			err = MPI_Isend(tmpbuf, sendcount, datatype,
				child[i], BLOCK_PARALLEL_BCAST, comm,
				&send_reqs[i]);
			if (err != MPI_SUCCESS) { line = __LINE__; goto error_hndl; }
		}

		err = MPI_Waitall(ct.nextsize, send_reqs,
			MPI_STATUSES_IGNORE);
		if (err != MPI_SUCCESS) { line = __LINE__; goto error_hndl; }
	}
	else {
		/*
		   Receive all segments from parent in a loop:
		   1) post irecv for the first segment
		   2) for segments 1 .. num_segments
		   - post irecv for the next segment
		   - wait on the previous segment to arrive
		   3) wait for the last segment
		*/
		req_index = 0;
		err = MPI_Irecv(tmpbuf, count_by_segment, datatype,
			ct.parent.rank, BLOCK_PARALLEL_BCAST, comm, &recv_reqs[req_index]);

		if (err != MPI_SUCCESS) { line = __LINE__; goto error_hndl; }

		for (segindex = 1; segindex < num_segments; segindex++) {
			req_index = req_index ^ 0x1;
			tmpbuf += realsegsize;
			/* post receive for the next segment */
			err = MPI_Irecv(tmpbuf + realsegsize, count_by_segment,
				datatype, ct.parent.rank, BLOCK_PARALLEL_BCAST,
				comm, &recv_reqs[req_index]);
			if (err != MPI_SUCCESS) { line = __LINE__; goto error_hndl; }
			/* wait on the previous segment */
			err = MPI_Wait(&recv_reqs[req_index ^ 0x1],
				MPI_STATUS_IGNORE);
			if (err != MPI_SUCCESS) { line = __LINE__; goto error_hndl; }
		}

		err = MPI_Wait(&recv_reqs[req_index], MPI_STATUS_IGNORE);
		if (err != MPI_SUCCESS) { line = __LINE__; goto error_hndl; }
	}

	return MPI_SUCCESS;

error_hndl:
	printf("%s:%4d\t Error occurred %d, rank %d\n",
		__FILE__, line, err, my_rank);
	(void)line;  // silence compiler warnings
	
	free_reqs(recv_reqs, 2);
	
	if (NULL != send_reqs) {
		free_reqs(send_reqs, ct.nextsize);
	}

	return err;
}