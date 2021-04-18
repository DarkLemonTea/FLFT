#include <mpi.h>
#include <stdio.h>

template <typename T>
struct AnyArray{
    private:
        T* Array;
        int Count;
    public:
        AnyArray(T*array, int count){
            Array = array;
            Count = count;
        }
		AnyArray(void**empty_array, int count){
			*empty_array = malloc(sizeof(T)*count);
			Array = (T*)(*empty_array);
			Count = count;
		}
		void set(T value){
			for(int n=0;n<Count;n++) Array[n] = value;
		}
        void MAX(T* Other){
            for(int n=0;n<Count;n++) if(Array[n] < Other[n]) Array[n]=Other[n];
        }
        void MIN(T* Other){
            for(int n=0;n<Count;n++) if(Array[n] > Other[n]) Array[n]=Other[n];
        }
        void SUM(T* Other){
            for(int n=0;n<Count;n++) Array[n]+=Other[n];
        }
        void PROD(T* Other){
            for(int n=0;n<Count;n++) Array[n]*=Other[n];
        }
        void BAND(T* Other){
            for(int n=0;n<Count;n++) Array[n]&=Other[n];
        }
        void BOR(T* Other){
            for(int n=0;n<Count;n++) Array[n]|=Other[n];
        }
        void BXOR(T* Other){
            for(int n=0;n<Count;n++) Array[n]^=Other[n];
        }
        void Op_reduce(MPI_Op op, void* other){
            switch(op){
                case MPI_MAX:  this->MAX((T*)other);  break;
                case MPI_MIN:  this->MIN((T*)other);  break;
                case MPI_SUM:  this->SUM((T*)other);  break;
                case MPI_PROD: this->PROD((T*)other); break;
                case MPI_BAND: this->BAND((T*)other); break;
                case MPI_BOR:  this->BOR((T*)other);  break;
                case MPI_BXOR: this->BXOR((T*)other); break;
                default: printf("\tReduce Error: Operator Unsupported.\n");
            }
        }
};

template<>
void AnyArray<float>::BAND(float* Other){printf("\tReduce Error: Float can not be BAND.\n");}
template<>
void AnyArray<float>::BOR(float* Other){printf("\tReduce Error: Float can not be BOR.\n");}
template<>
void AnyArray<float>::BXOR(float* Other){printf("\tReduce Error: Float can not be BXOR.\n");}
template<>
void AnyArray<double>::BAND(double* Other){printf("\tReduce Error: Double can not be BAND.\n");}
template<>
void AnyArray<double>::BOR(double* Other){printf("\tReduce Error: Double can not be BOR.\n");}
template<>
void AnyArray<double>::BXOR(double* Other){printf("\tReduce Error: Double can not be BXOR.\n");}

void local_reduce(MPI_Op op, void *source, void *target, int count, MPI_Datatype datatype){
    //º”µΩtarget¿Ô
	switch(datatype){
		case MPI_CHAR:   AnyArray<char>((char*)target,count).Op_reduce(op,source);   break;
		case MPI_SHORT:  AnyArray<short>((short*)target,count).Op_reduce(op,source);  break;
		case MPI_INT:    AnyArray<int>((int*)target,count).Op_reduce(op,source); break;
		case MPI_LONG:   AnyArray<long>((long*)target,count).Op_reduce(op,source);   break;
		case MPI_FLOAT:  AnyArray<float>((float*)target,count).Op_reduce(op,source);  break;
		case MPI_DOUBLE: AnyArray<double>((double*)target,count).Op_reduce(op,source); break;
		default: printf("\tReduce Error: Datatype Unsupported.\n");
	}
}