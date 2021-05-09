#ifndef HEADER_FILE
#define HEADER_FILE
#include"../fl_head.h"
#endif
#include"ring_connect.h"
#include"ring_multigather_lagging_procs.h"

//=============================================================================================================================
//lagging进程的主动阻塞状态，被捡回后需要确定一个rescuer
//=============================================================================================================================
void ring_block_and_rescuer_confirm(
	MPI_Comm comm,
	FD_var fd,
	Comm_proc *rescuer
) {
	int my_rank, comm_size;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &comm_size);
	printf("rank %d block!\n", my_rank);

	double block_time = 10 * fd.T_wait;
	Comm_proc savior, tmp_proc;
	savior.comm_stage = CONNECT_LAGGING;
	double cost_time;
	struct timeval st, start, time;

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

	init_ring(comm_size, my_rank, detector_stage, (*sp).Lagging_procs, &(*sp).ring);

	int result, res;
	result = ring_procs_connect(detector_stage, fd, comm, (*sp).Lagging_procs, &(*sp).ring);
	//printf("rank %d gets through the FD\n", my_rank);
	
	switch (result)
	{
	case RING_CONNECT_SUCCESS: {
		//直接通过
		res = FD_SUCCESS;
		break;
	}
	case RING_CONNECT_THROUGH: {
		//多聚集lagging进程
		multigather_lagging_procs(comm, (*sp).ring, &(*sp).Lagging_procs);
		//更新下邻居
		init_ring(comm_size, my_rank, detector_stage, (*sp).Lagging_procs, &(*sp).ring);
		res = FD_PASS;
		break;
	}
	case RING_CONNECT_LAGGING: {
		//printf("rank %d rescuer is %d\n", my_rank, rescuer.rank);
		ring_block_and_rescuer_confirm(comm, fd, &rescuer);
		
		if (rescuer.comm_stage == CONNECT_TIMEOUT_BREAK) {
			res = FD_FAILURE;
			break;
		}

		printf("rank %d waiting for rank %d to pick itself up\n", my_rank, rescuer.rank);
		while (1)
		{
			recv_procs_template(comm, &rescuer,
				RING_REVIVE_DATA_LAGGING_NUM, RING_REVIVE_DATA_LAGGING_PROCS,
				&(*sp).Lagging_procs);
			if (rescuer.comm_stage == FINISH) {
				(*sp).current_stage = rescuer.recv;
				break;
			}
		}

		init_ring(comm_size, my_rank, detector_stage, (*sp).Lagging_procs, &(*sp).ring);
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

