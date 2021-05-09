#pragma once

#include "../FLFT/Ring_FD/ring_FD.h"
#include "../FLFT/RT_FD/rt_FD.h"

//Failure detection
int FLFT_FD(
	//input
	int detector_stage,
	FD_var fd,
	MPI_Comm comm,
	//output
	Detector *sp //lagging_set
) {
	int res;
	switch ((*sp).type)
	{
	case TYPE_NONE: {
		res = FD_SUCCESS;
		break;
	}
	case TYPE_RING: {
		res = ring_FD(detector_stage, fd, comm, sp);
		break;
	}
	case TYPE_RING_TREE: {
		res = ring_tree_FD(detector_stage, fd, comm, sp);
		break;
	}
	default:
		res = ERROR;
		break;
	}

	return res;
}

