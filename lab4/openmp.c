#include "utility.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>
#include "numgen.c"

int checkIfPrime(unsigned long int inputArgument){
  if (inputArgument == 1 || inputArgument == 0)
    return 0;
  if (inputArgument < 4)
    return 1;
  long range = inputArgument / 2;
  for (long i = 2 ; i <= range;  i++){
    if (inputArgument % i == 0 ){
      return 0;
    }
  }
  return 1;
}


int main(int argc,char **argv) {


  Args ins__args;
  parseArgs(&ins__args, &argc, argv);

  //set number of threads
  omp_set_num_threads(ins__args.n_thr);
  
  //program input argument
  long inputArgument = ins__args.arg; 
  unsigned long int *numbers = (unsigned long int*)malloc(inputArgument * sizeof(unsigned long int));
  numgen(inputArgument, numbers);

  struct timeval ins__tstart, ins__tstop;
  gettimeofday(&ins__tstart, NULL);
  
  // run your computations here (including OpenMP stuff)
  unsigned long int result = 0;

  #pragma omp parallel
  {
    unsigned long int temp_result = 0;

    #pragma omp for
    for (unsigned long int i = 0; i < inputArgument; i++) {
      temp_result += checkIfPrime(numbers[i]);
    }

    #pragma omp critical
    {
      result += temp_result;
    }
  }
  
  printf("Liczb pierwszych: %ld\n", result);
  // synchronize/finalize your computations
  gettimeofday(&ins__tstop, NULL);
  ins__printtime(&ins__tstart, &ins__tstop, ins__args.marker);

}
