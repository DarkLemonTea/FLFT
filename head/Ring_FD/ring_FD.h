#include"ring_head.h"
#include"ring_connect.h"
#include"ring_gather_lagging_procs.h"
#include"ring_retrieve.h"

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
		if (clear.flag == 1 ) {
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

int ring_FD(
	//input
	int detector_stage,
	FD_var fd,
	MPI_Comm comm,
	//output
	Detector *sp //lagging_set
) {
	int my_rank, comm_size;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &comm_size);

	Comm_proc rescuer; //最终拯救者
	rescuer.rank = -1;

	init_ring(comm_size, my_rank, detector_stage, (*sp).Lagging_procs.num, (*sp).Lagging_procs.procs, &(*sp).ring);

	int result, res;
	result = ring_procs_connect(detector_stage, fd, comm, (*sp).Lagging_procs, &(*sp).ring);
	//printf("rank %d gets through the FD\n", my_rank);

	(*sp).Last_lagging_procs.num = (*sp).Lagging_procs.num;
	(*sp).Last_lagging_procs.procs = (*sp).Lagging_procs.procs;
	
	switch (result)
	{
	case RING_CONNECT_SUCCESS: {
		//直接通过
		res = FD_SUCCESS;
		break;
	}
	case RING_CONNECT_THROUGH: {
		//多聚集lagging进程
		multigather_lagging_procs(comm, (*sp).ring, &(*sp).Lagging_procs.num, &(*sp).Lagging_procs.procs);
		//更新下邻居
		init_ring(comm_size, my_rank, detector_stage, (*sp).Lagging_procs.num, (*sp).Lagging_procs.procs, &(*sp).ring);
		res = FD_SUCCESS;
		break;
	}
	case RING_CONNECT_LAGGING: {
		//printf("rank %d rescuer is %d\n", my_rank, rescuer.rank);
		block_and_rescuer_confirm(comm, fd, 
			(*sp).Lagging_procs.num, (*sp).Lagging_procs.procs, &rescuer);
		
		if (rescuer.comm_stage == CONNECT_TIMEOUT_BREAK) {
			res = FD_FAILURE;
			break;
		}

		printf("rank %d waiting for rank %d to pick itself up\n", my_rank, rescuer.rank);
		while (1)
		{
			recv_procs_template(comm, &rescuer,
				RING_REVIVE_DATA_LAGGING_NUM, RING_REVIVE_DATA_LAGGING_PROCS,
				&(*sp).Lagging_procs.num, &(*sp).Lagging_procs.procs);
			if (rescuer.comm_stage == FINISH) {
				(*sp).current_stage = rescuer.recv;
				break;
			}
		}

		init_ring(comm_size, my_rank, detector_stage, (*sp).Lagging_procs.num, (*sp).Lagging_procs.procs, &(*sp).ring);
		//清空消息缓存
		clear_message0(comm, RING_CONNECT1, (*sp).ring);
		clear_message2(comm, RING_CONNECT2);
		clear_message2(comm, RING_CONNECT3);
		clear_message1(comm, RING_RECONNECT1);
		clear_message2(comm, RING_RECONNECT2);
		clear_message2(comm, RING_RECONNECT3);
		clear_message1(comm, RING_LAGGING_PROBE);
		clear_message1(comm, RING_LAGGING_RESPOND);
		clear_message1(comm, RING_MISCONNECTION);
		//printf("recv stage is %d\n", (sp->current_stage));

		res = FD_REVIVE;
		break;
	}
	case RING_CONNECT_FAILURE: {
		res = FD_FAILURE;
		break;
	}
	default:
		break;
	}

	return res;
}

