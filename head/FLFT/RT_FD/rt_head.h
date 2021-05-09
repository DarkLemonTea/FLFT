#pragma once

#include"../Ring_FD/ring_head.h"

//=============================================================================================================================
//标识管理
//=============================================================================================================================

//返回结果
#define RT_CONNECT_SUCCESS 0
#define RT_CONNECT_THROUGH 1
#define RT_CONNECT_FAILURE 2
#define RT_CONNECT_LAGGING 3

//————————————————————————————————————————————————————————
//树方向
#define TREE_CONNECT1 10010
#define TREE_CONNECT2 10011
#define TREE_CONNECT3 10012
#define TREE_CONNECT4 10013

//计数
#define RING_COUNT_REDUCE 11000
#define RING_COUNT_BCAST 11001

#define TREE_COUNT_REDUCE 11002
#define TREE_COUNT_BCAST 11003

#define MISSING_REDUCE_COUNT_REQUEST 11004
#define MISSING_BCAST_COUNT_REQUEST 11005

//————————————————————————————————————————————————————————
//迟到进程处理
#define TREE_GATHER_LAGGING_NUM 20042
#define TREE_GATHER_LAGGING_PROCS 20043
#define TREE_MULTICAST_LAGGING_NUM 20044
#define TREE_MULTICAST_LAGGING_PROCS 20045
#define MISSING_GATHER_DATA_REQUEST 20046
#define MISSING_MULTICAST_DATA_REQUEST 20047

#define TREE_GATHER_REVIVE_NUM 20048
#define TREE_GATHER_REVIVE_PROCS 20049
#define TREE_MULTICAST_REVIVE_NUM 20050
#define TREE_MULTICAST_REVIVE_PROCS 20051
#define TREE_REVIVE_DATA_LAGGING_NUM 20052
#define TREE_REVIVE_DATA_LAGGING_PROCS 20053

//————————————————————————————————————————————————————————

#define TREE_ROOT 0
#define TREE_BRANCH 1
#define TREE_LEAF 2

//=============================================================================================================================
//Ring-Tree 等数据结构定义
//=============================================================================================================================

typedef struct tree {
	Comm_proc parent;
	Comm_proc left_child;
	Comm_proc right_child;
}Tree;

typedef struct last_connection_status {
	int p_lcs;
	int l_lcs;
	int r_lcs;
}Last_con_sta;

typedef struct ring_tree {
	int ring_scale;
	int ring_num;
	Ring ring;
	Tree tree;
	Last_con_sta lcs;
	int ide;       //根据树的位置，判断是根、枝、叶
}Ring_Tree;

//=============================================================================================================================
//Ring-Tree初始化操作
//=============================================================================================================================

void init_lcs(Last_con_sta *lcs) {
	(*lcs).p_lcs = 0;
	(*lcs).l_lcs = 0;
	(*lcs).r_lcs = 0;
}


//寻找环内最近有效对象
int left_rt_proc(int comm_size, int my_rank, int ring_scale) {
	int order = my_rank / ring_scale; //环规模
	int ring_rank = my_rank % ring_scale;  //进程在环内的序号
	int ring_num;  //环进程数量
	if (order == comm_size / ring_scale) {
		ring_num = comm_size % ring_scale;
	}
	else { ring_num = ring_scale; }
	int left_rt = (ring_num + ring_rank - 1) % ring_num;
	return left_rt + ring_scale * order;
}

int right_rt_proc(int comm_size, int my_rank, int ring_scale) {
	int order = my_rank / ring_scale; //环规模
	int ring_rank = my_rank % ring_scale;  //进程在环内的序号
	int ring_num;  //环进程数量
	if (order == comm_size / ring_scale) {
		ring_num = comm_size % ring_scale;
	}
	else { ring_num = ring_scale; }
	int right_rt = (ring_rank + 1) % ring_num;
	return right_rt + ring_scale * order;
}

int find_rt_valid(
	int comm_size,
	int my_rank,
	int ring_scale,
	P_set rlp,
	char l_or_r
) {
	int tmp_proc;
	if (l_or_r == 'l') {
		tmp_proc = left_rt_proc(comm_size, my_rank, ring_scale);
		while (in(rlp, tmp_proc)) {
			tmp_proc = left_rt_proc(comm_size, tmp_proc, ring_scale);
		}
	}
	else if (l_or_r == 'r') {
		tmp_proc = right_rt_proc(comm_size, my_rank, ring_scale);
		while (in(rlp, tmp_proc)) {
			tmp_proc = right_rt_proc(comm_size, tmp_proc, ring_scale);
		}
	}
	return tmp_proc;
}

int rt_parent(int comm_size, int my_rank, int ring_scale) {
	int parent;
	int tree_num = my_rank % ring_scale; //树编号，第几棵树
	int tree_rank = my_rank / ring_scale; //树中所处的位置
	if (tree_rank == 0) {
		parent = EMPTY;
	}
	else {
		parent = (tree_rank / 2) - (1 - tree_rank % 2);
		parent = parent * ring_scale + tree_num;
	}
	return parent;
}

int rt_lchild(int comm_size, int my_rank, int ring_scale) {
	int lchild;
	int tree_num = my_rank % ring_scale; //树编号，第几棵树
	int tree_rank = my_rank / ring_scale; //树中所处的位置

	lchild = (tree_rank * 2 + 1)*ring_scale + tree_num;
	if (lchild >= comm_size) { lchild = EMPTY; }
	return lchild;
}

int rt_rchild(int comm_size, int my_rank, int ring_scale) {
	int rchild;
	int tree_num = my_rank % ring_scale; //树编号，第几棵树
	int tree_rank = my_rank / ring_scale; //树中所处的位置

	rchild = (tree_rank * 2 + 2)*ring_scale + tree_num;
	if (rchild >= comm_size) { rchild = EMPTY; }
	return rchild;
}

void init_ring_tree(
	int comm_size,
	int my_rank,
	int ring_scale,
	int detector_stage,
	P_set rlp,
	Ring_Tree *rt
) {
	(*rt).ring_scale = ring_scale;

	if (my_rank / ring_scale == comm_size / ring_scale) {
		(*rt).ring_num = comm_size % ring_scale;
	}
	else { (*rt).ring_num = ring_scale; }

	(*rt).ring.left_proc.rank = find_rt_valid(comm_size, my_rank, ring_scale, rlp, 'l');
	(*rt).ring.right_proc.rank = find_rt_valid(comm_size, my_rank, ring_scale, rlp, 'r');

	(*rt).tree.parent.rank = rt_parent(comm_size, my_rank, ring_scale);
	(*rt).tree.left_child.rank = rt_lchild(comm_size, my_rank, ring_scale);
	(*rt).tree.right_child.rank = rt_rchild(comm_size, my_rank, ring_scale);

	(*rt).ring.left_proc.comm_stage = CONNECT_STATE1;
	(*rt).ring.right_proc.comm_stage = CONNECT_STATE1;

	if ((*rt).tree.parent.rank == EMPTY) { (*rt).tree.parent.comm_stage = CONNECT_FINISH; }
	else { (*rt).tree.parent.comm_stage = CONNECT_STATE1; }

	if ((*rt).tree.left_child.rank == EMPTY) { (*rt).tree.left_child.comm_stage = CONNECT_FINISH; }
	else { (*rt).tree.left_child.comm_stage = CONNECT_STATE1; }

	if ((*rt).tree.right_child.rank == EMPTY) { (*rt).tree.right_child.comm_stage = CONNECT_FINISH; }
	else { (*rt).tree.right_child.comm_stage = CONNECT_STATE1; }

	(*rt).ring.left_proc.send = detector_stage;
	(*rt).ring.right_proc.send = detector_stage;
	(*rt).tree.parent.send = detector_stage;
	(*rt).tree.left_child.send = detector_stage;
	(*rt).tree.right_child.send = detector_stage;

	//init_lcs(&(*rt).lcs);

	//printf("rank %d, ring_scale is %d, ring l %d,ring r is %d\t tree p %d, tree l %d, tree r %d\n",
	//	my_rank, ring_scale, (*rt).ring.left_proc.rank, (*rt).ring.right_proc.rank,
	//	(*rt).tree.parent.rank, (*rt).tree.left_child.rank, (*rt).tree.right_child.rank);
}

//统计迟到进程（进程自身负责的部分）
void rt_statistics_lagged_procs(
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