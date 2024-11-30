#include "utility.h"
#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <sys/time.h>
#include "numgen.c"

__host__
void errorexit(const char *s) {
    printf("\n%s",s);	
    exit(EXIT_FAILURE);	 	
}

__global__ 
void checkIfPrime(unsigned long int* inputArgument, unsigned long int* returnArg){
  unsigned long int index=blockIdx.x*blockDim.x+threadIdx.x;
  unsigned long int my_number = (unsigned long int)inputArgument[index];
  if (my_number == 1 || my_number == 0) {
    returnArg[index] = 0;
    return;
  }

  if (my_number < 4){
    returnArg[index] = 1;
    return;
  }
    
  long range = my_number / 2 + 1;
  for (long i = 2 ; i*i <= my_number;  i++){
    if (my_number % i == 0 ){
      returnArg[index] = 0;
      return;
    }
  }
  returnArg[index] = 1;
  return;
}



int main(int argc,char **argv) {

  Args ins__args;
  parseArgs(&ins__args, &argc, argv);

  int threadsinblock=1024;
  int blocksingrid=10000;	
  long long result;
  long long size = threadsinblock*blocksingrid;
  
  //program input argument
  long inputArgument = ins__args.arg; 
  unsigned long int *numbers = (unsigned long int*)malloc(size * sizeof(unsigned long int));
  for (int i =0;i<size;i++)
    numbers[i]=0;

  numgen(inputArgument, numbers);

  struct timeval ins__tstart, ins__tstop;
  gettimeofday(&ins__tstart, NULL);
  
  // run your CUDA kernel(s) here

    
    
    //memory allocation on host
    unsigned long int *hresults=(unsigned long int*)malloc(size*sizeof(unsigned long int));
    if (!hresults) errorexit("Error allocating memory on the host");	

    unsigned long int *dresults=NULL;
    unsigned long int *dnumbers=NULL;

    if (cudaSuccess!=cudaMalloc((void **)&dresults,size*sizeof(unsigned long int)))
      errorexit("Error allocating memory on the GPU");

    if (cudaSuccess!=cudaMalloc((void **)&dnumbers,size*sizeof(unsigned long int)))
      errorexit("Error allocating memory on the GPU");

    if (cudaSuccess!=cudaMemcpy(dnumbers,numbers,size*sizeof(unsigned long int),cudaMemcpyHostToDevice))
       errorexit("Error copying numbers");

    checkIfPrime<<<blocksingrid,threadsinblock>>>(dnumbers, dresults);
    if (cudaSuccess!=cudaGetLastError())
      errorexit("Error during kernel launch");
  
    if (cudaSuccess!=cudaMemcpy(hresults,dresults,size*sizeof(unsigned long int),cudaMemcpyDeviceToHost))
       errorexit("Error copying results ");


    //calculate sum of all elements on CPU side
    result=0;

    for(int i=0;i<size;i++) {
      result += hresults[i];
    }

    printf("\nThe final result is %lld\n",result);

    // synchronize/finalize your CUDA computations

  gettimeofday(&ins__tstop, NULL);
  ins__printtime(&ins__tstart, &ins__tstop, ins__args.marker);


}
