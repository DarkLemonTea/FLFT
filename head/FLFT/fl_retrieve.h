#pragma once


#include"RT_FD/rt_retrieve.h"

void retrieve_procs(
	int detector_stage,
	MPI_Comm comm,
	double T_retrieve,
	Detector *sp
) {
	switch ((*sp).type)
	{
	case TYPE_NONE: {
		break;
	}
	case TYPE_RING: {
		ring_retrieve_procs(detector_stage, comm, T_retrieve, sp);
		break;
	}
	case TYPE_RING_TREE: {
		rt_retrieve_procs(detector_stage, comm, T_retrieve, sp);
		break;
	}
	default:
		break;
	}
}