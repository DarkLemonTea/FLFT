# Ring Failure Detection

We use function Ring_FD to implement the failure detection.

```
int ring_FD( int detector_stage, FD_var fd, MPI_Comm comm, Detector *sp)
```

- The **detector stage** is an increasing sequence number, which can be used to judge the failure of the process and the reuse after the lagging processes is retrieved.

- **FD_var** is the a variable used to set parameters, including wall time, wait time, release rate, and more.

- **Detector** is used to store pointers to sets and variables, including Lagging processes sets, Retrieved processes sets and more.

For detailed usage, please refer to the case in ../src/test.c
