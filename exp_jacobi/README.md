# Block jacobi iterations

Block jacobi with FT compares with Block jacobi without FT

The result of experiments on Tianhe2:

Each process is responsible for the calculation of a 3000*3000 matrix block

We use the ring failure detection.

As the number of processes increases, the efficiency will be reduced.

Next we will use Ring-Tree failure detection to improve efficiency.


The iteration converges at step 784

- **bj-nfd:** use time of block-jacobi pure calculation
- **bj-fd-nf:** use time of block-jacobi with failure detection, backup and fault tolerance
- **bj-fd-sf:** use time of fault-tolerant block jacobi in the case of simulated failure.

*(The unit of time is milliseconds.)*

node | processes | bj-nfd | bj-fd-nf | bj-fd-sf
---|---|---|---|---
6 | 144 | 17141.000 | 18960.530  | 21655.565
24 | 576 | 17280.531 | 	21051.824 | 24010.787
43 | 1024 | 18538.327000 | 22777.632000 | 	26948.070000
96 | 2304 | 19680.548000 | 27453.336000 | 	31797.492000
150 | 3600 | 19211.530000 | 34977.161000 | 	36780.654000

