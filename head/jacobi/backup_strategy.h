#include<stdio.h>
#include<stdlib.h>

#ifndef HEADER_FILE
#define HEADER_FILE
#include "../Ring_FD/ring_head.h"
#include "block_parallel.h"
#endif

/*
���ڷֿ鲢�еı��ݴ洢��ָ�����
�ֿ鲢��jacobi�ص�
row reduce
col bcast
root����ʧ����ͬ�н��̳е�
���root�������ݱ�����ͬ�н���

��root����ʧ��
ͬ�н��̸���local reduce�������ͬ�н��̳��б���

root���̲����б���

���ݺͻ�ȡ��Ϊ������ͬ�ĺ�����ÿ�����±��ݺ󣬻Ḳ����һ�ε���Ч״̬
*/

#define REBUILD_SUCCESS 0
#define REBUILD_FAILURE 1

typedef struct backup_obj {
	int rank;
	Comm_proc des;
	Comm_proc src;
}Backup;

void back_up_index(
	int rank,
	int n,
	Backup *b
) {
	if (rank == EMPTY) {
		(*b).des.rank = EMPTY; 
		(*b).src.rank = EMPTY;
		return;
	}

	(*b).rank = rank;
	int des, src;

	//root
	if (rank / n == rank % n) {
		//root����ֻ��Ҫ�ҵ����ݽ��̣�����Ҫ�����������̵�����
		src = EMPTY;
		des = rank % n + ((rank / n + 1) % n)*n;
	
		init_proc(src, &(*b).src);
		init_proc(des, &(*b).des);
	}
	else {
		src = (rank%n - 1 + n) % n + (rank / n)*n;
		des = (rank%n + 1) % n + (rank / n)*n;

		while (src % n == src / n)
		{
			src = (src%n - 1 + n) % n + (src / n)*n;
		}
		while (des % n == des / n)
		{
			des = (des%n + 1) % n + (des / n)*n;
		}
		init_proc(des, &(*b).des);
		init_proc(src, &(*b).src);
	}
}

int rebuild_comm_graph(
	Detector fd,
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

	int i = 0, index, num = n, ll_num = 0, ind = 0, tmp_rank;
	int *procs, *ll_procs;

	Backup tmp;
	
	if (row_or_col == 'r') {
		//��rΪ��ʱ���ݴ���ͨ����reduce��
		ind = 0;
		//���ȴ�lagging������ɸѡ��ͬ�е�lagging����
		ll_procs = (int*)calloc(fd.Lagging_procs.num, sizeof(int));
		for (i = 0; i < fd.Lagging_procs.num; i++) {
			if (fd.Lagging_procs.procs[i] / n == row_ind) {
				ll_procs[ind] = fd.Lagging_procs.procs[i];
				ind += 1;
			}
		}
		ll_num = ind;
		num = n - ll_num;

		procs = (int*)calloc(num, sizeof(int));

		procs[0] = root;
		if (in(ll_num, ll_procs, root)) {
			back_up_index(root, n, &tmp);
			if (in(fd.Lagging_procs.num, fd.Lagging_procs.procs, tmp.des.rank)) {
				if (ll_procs != NULL) { free(ll_procs); }
				if (procs != NULL) { free(procs); }
				return REBUILD_FAILURE;
			}
			else {
				//�ҵ��������
				procs[0] == tmp.des.rank;
			}
		}

		ind = 1;

		for (i = 1; i < n; i++) {
			tmp_rank = row_ind * n + i;
			if (in(ll_num, ll_procs, tmp_rank)) {
				back_up_index(tmp_rank, n, &tmp);
				//����޷��ָ������ж��˳�
				if (in(fd.Lagging_procs.num, fd.Lagging_procs.procs, tmp.des.rank)) {
					if (ll_procs != NULL) { free(ll_procs); }
					if (procs != NULL) { free(procs); }
					return REBUILD_FAILURE;
				}
				else {
					continue;
				}		
			}
			else {
				procs[ind] = tmp_rank;
				if (tmp_rank == root) {
					procs[ind] = row_ind * n;
				}
				if (my_rank == procs[ind]) {
					index = ind;
				}
				ind += 1;
			}		 
		}

		if (my_rank == root) { index = 0; }
	}
	else if (row_or_col = 'c') {
		//��rΪ��ʱ���ݴ���ͨ����reduce��
		ind = 0;
		//���ȴ�lagging������ɸѡ��ͬ�е�lagging����
		ll_procs = (int*)calloc(fd.Lagging_procs.num, sizeof(int));
		for (i = 0; i < fd.Lagging_procs.num; i++) {
			if (fd.Lagging_procs.procs[i] % n == col_ind) {
				ll_procs[ind] = fd.Lagging_procs.procs[i];
				ind += 1;
			}
		}
		ll_num = ind;
		num = n;

		procs = (int*)calloc(num, sizeof(int));

		procs[0] = root;
		if (in(ll_num, ll_procs, root)) {
			back_up_index(root, n, &tmp);
			if (in(fd.Lagging_procs.num, fd.Lagging_procs.procs, tmp.des.rank)) {
				if (ll_procs != NULL) { free(ll_procs); }
				if (procs != NULL) { free(procs); }
				return REBUILD_FAILURE;
			}
			else {
				//�ҵ��������
				procs[0] == tmp.des.rank;
			}
		}

		ind = 1;
		
		for (i = 1; i < num; i++) {
			tmp_rank = col_ind + i * n;
			if (in(ll_num, ll_procs, tmp_rank)) {
				back_up_index(tmp_rank, n, &tmp);
				//����޷��ָ������ж��˳�
				if (in(fd.Lagging_procs.num, fd.Lagging_procs.procs, tmp.des.rank)) {
					if (ll_procs != NULL) { free(ll_procs); }
					if (procs != NULL) { free(procs); }
					return REBUILD_FAILURE;
				}
				else {
					procs[ind] = tmp.des.rank;
					ind += 1;
				}
			}
			else {
				procs[ind] = tmp_rank;
				if (my_rank == tmp_rank) {
					index = ind;
				}
				if (tmp_rank%n == tmp_rank/n) {
					procs[ind] = col_ind;
					if (my_rank == col_ind) {
						index = ind;
					}
					if (in(fd.Lagging_procs.num, fd.Lagging_procs.procs, col_ind)) {
						back_up_index(col_ind, n, &tmp);
						procs[ind] = tmp.des.rank;
					}
				}
				ind += 1;
			}
		}

		if (my_rank == root) { index = 0; }
	}

	//�������꣬��ȡ��Ӧ����ֵ
		//������
	if (index == 0) { (*ct).parent.rank = -1; }
	else { (*ct).parent.rank = procs[(index - 1) / 2]; }
	//����
	if (index * 2 + 1 >= num) { (*ct).lchild.rank = -1; }
	else { (*ct).lchild.rank = procs[index * 2 + 1]; }
	//�Һ���
	if (index * 2 + 2 >= num) { (*ct).rchild.rank = -1; }
	else { (*ct).rchild.rank = procs[index * 2 + 2]; }

	if (ll_procs != NULL) { free(ll_procs); }
	if (procs != NULL) { free(procs); }

	if ((*ct).lchild.rank > -1) {
		(*ct).nextsize += 1;
	}
	if ((*ct).rchild.rank > -1) {
		(*ct).nextsize += 1;
	}

	return REBUILD_SUCCESS;
}
