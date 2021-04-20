# FLFT
User-level failure detection and fault tolerant algorithm based on fail-lagging model

The Fail-lagging model is a failure model under the loose BSP model. We can treat each detection as a loose barrier, which does not require all processes to arrive before the arrival processes can be released, and a detection mode independent of the application algorithm is implemented at the application level.

Under the parallel system model, it is impossible to distinguish whether a failed process is crashing or very slow, the traditional implementation method based on heartbeat detection is through all-to-all heartbeat monitoring. When the number of processes is large, this method will cause great communication overhead. At the same time, it will reduce the operating efficiency of the application. It is difficult to select the timeout value of the Gossip method, and there is no way to avoid the false negative problem. Based on ring observation and hypercube failures propagation, false negatives or continuous processes failed cannot be dealt with.

In fact, it is the contradiction between local information and global information. To obtain accurate results, huge communication overhead is required. If you want to reduce communication overhead, it is difficult to grasp the overall information. This is a dilemma.

False negative processes often cause the application to run ineffectively, with no response and no output for a long time, without reporting an error or interrupting the exit, which causes a serious waste of calculations. In order to avoid this problem during the running of the application, we propose to use the Fail-lagging model to describe the failures. 

Under the loose BSP model, after the number of arrival processes reaches a certain number that arrow arrival processes to pass and waits for enough time, it is no longer necessary to wait for the non-arrival process, and pass the barrier directly with fault-tolerant to execute the next task.

Failure detection uses timing and counting methods, local connections establish connected maximally connected subgraph, count on the maximally connected subgraph, and use multi-reduce communication to quickly obtain global information, without the need to send all to all heartbeat, which greatly reduces the complexity of collective communication and can effectively avoid false negative results.

However, with this judgment mode, the timeout value is usually not extended too much, and there may be available processes that become isolated points, which is a false positive. The Fail-lagging model does not immediately terminate the lagging process. If the lagging process is still available, it can be retrieved and reused.

The implementation of failure detection is based on ring, ring-tree, ring-butterfly. Collective communication uses point-to-point communication to improve the realization, increase the failure tolerance function, can reconstruct the communication structure according to the lagging process set, and ensure that the collective communication can be completed.

The experiments and results can be found in the folder prefixed with "exp_".

The experimental results of the Jacobi experiment have been displayed in the folder "exp_jacobi"