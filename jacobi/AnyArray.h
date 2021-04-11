#include<stdio.h>

typedef struct AnyArray{
    //数据域
    int Count;
    MPI_Datatype Dtype;
    char*   Array00;
    short*  Array01;
    int*    Array02;
    long*   Array03;
    float*  Array09;
    double* Array10;
    //方法
    AnyArray(void*source, MPI_Datatype datatype, int count){
        Dtype = datatype;
        Count = count;
        switch(datatype){
            case MPI_CHAR:   Array00 = (char*)  source; break;
            case MPI_SHORT:  Array01 = (short*) source; break;
            case MPI_INT:    Array02 = (int*)   source; break;
            case MPI_LONG:   Array03 = (long*)  source; break;
            case MPI_FLOAT:  Array09 = (float*) source; break;
            case MPI_DOUBLE: Array10 = (double*)source; break;
        }
    }

    void MAX(AnyArray& other){
        switch(Dtype){
            case MPI_CHAR:   for(int n=0;n<Count;n++) if(Array00[n] < other.Array00[n]) Array00[n]=other.Array00[n]; break;
            case MPI_SHORT:  for(int n=0;n<Count;n++) if(Array01[n] < other.Array01[n]) Array01[n]=other.Array01[n]; break;
            case MPI_INT:    for(int n=0;n<Count;n++) if(Array02[n] < other.Array02[n]) Array02[n]=other.Array02[n]; break;
            case MPI_LONG:   for(int n=0;n<Count;n++) if(Array03[n] < other.Array03[n]) Array03[n]=other.Array03[n]; break;
            case MPI_FLOAT:  for(int n=0;n<Count;n++) if(Array09[n] < other.Array09[n]) Array09[n]=other.Array09[n]; break;
            case MPI_DOUBLE: for(int n=0;n<Count;n++) if(Array10[n] < other.Array10[n]) Array10[n]=other.Array10[n]; break;
        }
    }
    void MIN(AnyArray& other){
        switch(Dtype){
            case MPI_CHAR:   for(int n=0;n<Count;n++) if(Array00[n] > other.Array00[n]) Array00[n]=other.Array00[n]; break;
            case MPI_SHORT:  for(int n=0;n<Count;n++) if(Array01[n] > other.Array01[n]) Array01[n]=other.Array01[n]; break;
            case MPI_INT:    for(int n=0;n<Count;n++) if(Array02[n] > other.Array02[n]) Array02[n]=other.Array02[n]; break;
            case MPI_LONG:   for(int n=0;n<Count;n++) if(Array03[n] > other.Array03[n]) Array03[n]=other.Array03[n]; break;
            case MPI_FLOAT:  for(int n=0;n<Count;n++) if(Array09[n] > other.Array09[n]) Array09[n]=other.Array09[n]; break;
            case MPI_DOUBLE: for(int n=0;n<Count;n++) if(Array10[n] > other.Array10[n]) Array10[n]=other.Array10[n]; break;
        }
    }
    void SUM(AnyArray& other){
        switch(Dtype){
            case MPI_CHAR:   for(int n=0;n<Count;n++) Array00[n]+=other.Array00[n]; break;
            case MPI_SHORT:  for(int n=0;n<Count;n++) Array01[n]+=other.Array01[n]; break;
            case MPI_INT:    for(int n=0;n<Count;n++) Array02[n]+=other.Array02[n]; break;
            case MPI_LONG:   for(int n=0;n<Count;n++) Array03[n]+=other.Array03[n]; break;
            case MPI_FLOAT:  for(int n=0;n<Count;n++) Array09[n]+=other.Array09[n]; break;
            case MPI_DOUBLE: for(int n=0;n<Count;n++) Array10[n]+=other.Array10[n]; break;
        }
    }
    void PROD(AnyArray& other){
        switch(Dtype){
            case MPI_CHAR:   for(int n=0;n<Count;n++) Array00[n]*=other.Array00[n]; break;
            case MPI_SHORT:  for(int n=0;n<Count;n++) Array01[n]*=other.Array01[n]; break;
            case MPI_INT:    for(int n=0;n<Count;n++) Array02[n]*=other.Array02[n]; break;
            case MPI_LONG:   for(int n=0;n<Count;n++) Array03[n]*=other.Array03[n]; break;
            case MPI_FLOAT:  for(int n=0;n<Count;n++) Array09[n]*=other.Array09[n]; break;
            case MPI_DOUBLE: for(int n=0;n<Count;n++) Array10[n]*=other.Array10[n]; break;
        }
    }
    void BAND(AnyArray& other){
        switch(Dtype){
            case MPI_CHAR:   for(int n=0;n<Count;n++) Array00[n]&=other.Array00[n]; break;
            case MPI_SHORT:  for(int n=0;n<Count;n++) Array01[n]&=other.Array01[n]; break;
            case MPI_INT:    for(int n=0;n<Count;n++) Array02[n]&=other.Array02[n]; break;
            case MPI_LONG:   for(int n=0;n<Count;n++) Array03[n]&=other.Array03[n]; break;
            case MPI_FLOAT:  /* 不支持 */ break;
            case MPI_DOUBLE: /* 不支持 */ break;
        }
    }
    void BOR(AnyArray& other){
        switch(Dtype){
            case MPI_CHAR:   for(int n=0;n<Count;n++) Array00[n]|=other.Array00[n]; break;
            case MPI_SHORT:  for(int n=0;n<Count;n++) Array01[n]|=other.Array01[n]; break;
            case MPI_INT:    for(int n=0;n<Count;n++) Array02[n]|=other.Array02[n]; break;
            case MPI_LONG:   for(int n=0;n<Count;n++) Array03[n]|=other.Array03[n]; break;
            case MPI_FLOAT:  /* 不支持 */ break;
            case MPI_DOUBLE: /* 不支持 */ break;
        }
    }
    void BXOR(AnyArray& other){
        switch(Dtype){
            case MPI_CHAR:   for(int n=0;n<Count;n++) Array00[n]^=other.Array00[n]; break;
            case MPI_SHORT:  for(int n=0;n<Count;n++) Array01[n]^=other.Array01[n]; break;
            case MPI_INT:    for(int n=0;n<Count;n++) Array02[n]^=other.Array02[n]; break;
            case MPI_LONG:   for(int n=0;n<Count;n++) Array03[n]^=other.Array03[n]; break;
            case MPI_FLOAT:  /* 不支持 */ break;
            case MPI_DOUBLE: /* 不支持 */ break;
        }
    }
}AnyArray;
