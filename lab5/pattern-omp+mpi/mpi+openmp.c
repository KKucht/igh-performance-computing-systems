#include "utility.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>
#include <omp.h>
#include "numgen.c"
#include <stdbool.h>


#define DATA 10
#define RESULT 11
#define FINISH 12
#define ENDING 13
#define SIZE_OF_BLOCK 1000
#define FOR_RECIVE 0
#define FOR_SEND 1

unsigned long int next_block(unsigned long int x){
  return (x + 1) % 3;
}

unsigned long int checkHowMany( long int size,  long int nr){
  long int size_to_prepare = (nr+1)*SIZE_OF_BLOCK;
  if(size - size_to_prepare > 0 ){
    size_to_prepare = SIZE_OF_BLOCK;
  }
  else{
    size_to_prepare = size - (nr)*SIZE_OF_BLOCK;
  }
  return size_to_prepare;
}

void prepareBlock(unsigned long int* input, unsigned long int* output, unsigned long int size) {
  for (unsigned long int i = 0; i < SIZE_OF_BLOCK; i++) {
    if (i >= size) {
      output[i] = 0;
    }
    else{
      output[i] = input[i];
    }
  }
  
}

void firstPrepareBlock(unsigned long int* input, unsigned long int* output, unsigned long int size) {
  for (unsigned long int i = 0; i < 2* SIZE_OF_BLOCK; i++) {
    if (i >= size) {
      output[i] = 0;
    }
    else{
      output[i] = input[i];
    }
  }
  
  
}

int checkIfPrime(unsigned long int inputArgument){
  if (inputArgument == 1 || inputArgument == 0)
    return 0;
  if (inputArgument < 4)
    return 1;
  long range = inputArgument / 2 + 1;
  for (long i = 2 ; i*i <= inputArgument;  i++){
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

  struct timeval ins__tstart, ins__tstop;

  int threadsupport;
  int myrank,nproc;
  unsigned long int *numbers;
  // Initialize MPI with desired support for multithreading -- state your desired support level

  MPI_Init_thread(&argc, &argv,MPI_THREAD_FUNNELED,&threadsupport); 

  if (threadsupport<MPI_THREAD_FUNNELED) {
    printf("\nThe implementation does not support MPI_THREAD_FUNNELED, it supports level %d\n",threadsupport);
    MPI_Finalize();
    return -1;
  }
  
  // obtain my rank
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  // and the number of processes
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);

  MPI_Request *requests;
  unsigned long int* my_block;

  if(!myrank){
      	gettimeofday(&ins__tstart, NULL);
    numbers = (unsigned long int*)malloc(inputArgument * sizeof(unsigned long int));
  	numgen(inputArgument, numbers);

    int nr_of_blocks = nproc - 1;
    if (nr_of_blocks < 2)
      nr_of_blocks = 2;
    my_block = (unsigned long int*)malloc(nr_of_blocks*SIZE_OF_BLOCK * sizeof(unsigned long int));

    requests = (MPI_Request *) malloc (2 * (nproc - 1) * sizeof (MPI_Request));

    for (unsigned long int i = 0; i < 2 * (nproc - 1); i++)
  		requests[i] = MPI_REQUEST_NULL;
  }
  else{

    my_block = (unsigned long int*)malloc(3*SIZE_OF_BLOCK * sizeof(unsigned long int));
    //printf("requests %ld",2 *  sizeof (MPI_Request));
    requests = (MPI_Request *) malloc (2 *  sizeof (MPI_Request));

    for (unsigned long int i = 0; i < 2 ; i++)
  		requests[i] = MPI_REQUEST_NULL;
  }
  // run your computations here (including MPI communication and OpenMP stuff)

    MPI_Status status;
  unsigned long int my_result = 0;
  unsigned long int *my_flags_to_not_send;
  unsigned long int *result_temp;

  if(myrank == 0 ){
    int index = 0;
    int licznik = 0;
    int requestcompleted = -1; 
    unsigned long int size_to_prepare = 0;
    unsigned long int nr_of_pack = inputArgument / SIZE_OF_BLOCK;

    if (inputArgument % SIZE_OF_BLOCK != 0){
      nr_of_pack++;
    }
    
    my_flags_to_not_send = (unsigned long int *) malloc ((nproc - 1) *  sizeof (unsigned long int));
    for (unsigned long int i = 0; i < (nproc - 1); i++)
  		my_flags_to_not_send[i] = 0;;

    
    result_temp = (unsigned long int *) malloc ((nproc - 1) *  sizeof (unsigned long int));

    unsigned long int nr = 0;
    unsigned long int slave_proc = 1;
    unsigned long int nr_results = 0;

    while (slave_proc < nproc ) 
    {
      // send first 2 package
      size_to_prepare = checkHowMany(inputArgument, nr);
      nr++;
      size_to_prepare += checkHowMany(inputArgument, nr);
      // printf("%ld\n", size_to_prepare);
      firstPrepareBlock(&(numbers[(slave_proc - 1)*2*SIZE_OF_BLOCK]), my_block, size_to_prepare);
      MPI_Send(my_block,SIZE_OF_BLOCK*2,MPI_UNSIGNED_LONG,slave_proc,DATA,MPI_COMM_WORLD);
      nr++;
      // printf("Send first package to %ld\n",slave_proc);
      // fflush(stdout);
      slave_proc++;
    }
    // send first nonblocking
    for(unsigned long int i=1;i<nproc;i++) {
      // printf("Start Irecv for %ld\n", i);
      // fflush(stdout);
      MPI_Irecv(&(result_temp[ i - 1 ]),1,MPI_UNSIGNED_LONG,i,RESULT,MPI_COMM_WORLD,&(requests[i - 1]));
    }

    for(unsigned long int i=1;i<nproc;i++) {
      size_to_prepare = checkHowMany(inputArgument, nr);
      prepareBlock(&(numbers[nr*SIZE_OF_BLOCK]), &(my_block[(i - 1 )*SIZE_OF_BLOCK]), size_to_prepare);
      // printf("Start Isend for %ld\n", i);
      // fflush(stdout);
      MPI_Isend(&(my_block[(i - 1 ) * SIZE_OF_BLOCK]),SIZE_OF_BLOCK,MPI_UNSIGNED_LONG,i,DATA,MPI_COMM_WORLD, &(requests[(i - 1) + (nproc - 1)]));
      nr++;
    }

    do {
      // printf("master: czekam na jakiś recive\n");
      // fflush(stdout);
      MPI_Waitany ((nproc - 1), &(requests[0]), &requestcompleted, MPI_STATUS_IGNORE);
      // printf("master: dostałem revice od %d\n", requestcompleted);
      // printf("master: czekam na send od %d\n", requestcompleted);
      // fflush(stdout);
      MPI_Wait(&(requests[requestcompleted + nproc - 1]), MPI_STATUS_IGNORE);
      
      if (nr < nr_of_pack) {
        // printf("master: wysyłam mu paczkę %d\n", requestcompleted);
        // fflush(stdout);
        size_to_prepare = checkHowMany(inputArgument, nr);
        prepareBlock(&(numbers[nr*SIZE_OF_BLOCK]), &(my_block[requestcompleted*SIZE_OF_BLOCK]), size_to_prepare);
        MPI_Isend(&(my_block[requestcompleted*SIZE_OF_BLOCK]),SIZE_OF_BLOCK,MPI_UNSIGNED_LONG,requestcompleted + 1,DATA,MPI_COMM_WORLD,&(requests[requestcompleted+nproc - 1]));
        nr++;
      }
      else if (my_flags_to_not_send[requestcompleted] == 0){
        // printf("master: informuję o zakończeniu %d\n", requestcompleted);
        // fflush(stdout);
        my_flags_to_not_send[requestcompleted] = 1;
        MPI_Isend(&(my_block[requestcompleted*SIZE_OF_BLOCK]),SIZE_OF_BLOCK,MPI_UNSIGNED_LONG,requestcompleted + 1,FINISH,MPI_COMM_WORLD, &(requests[requestcompleted+nproc - 1]));
      }
      // printf("master: zbieram wynik %d\n", requestcompleted);
      // fflush(stdout);
      my_result += result_temp[requestcompleted];
      nr_results++;
      if (nr_results < nr_of_pack){
        // printf("master: znowu pytam %d\n", requestcompleted);
        // fflush(stdout);
        MPI_Irecv(&(result_temp[requestcompleted]),1,MPI_UNSIGNED_LONG,requestcompleted+1,RESULT,MPI_COMM_WORLD,&(requests[requestcompleted]));
      }
      else{
        for (unsigned long int i=1;i<nproc;i++){
          if (requestcompleted != i - 1){
            MPI_Cancel (&(requests[i - 1]));
          }
          
        }
        break;
      }

    } while (true);
    
    
    // printf ("master: wynik: %ld\n", my_result);
    // printf("master: koniec\n");
    // fflush(stdout);
    for (unsigned long int x = 1 ; x < nproc ;x++)
      MPI_Recv(result_temp,1,MPI_UNSIGNED_LONG,x,ENDING,MPI_COMM_WORLD,&status);
  }
  else {
    // printf("slave: start\n");
    bool my_recive_tag = false;
    unsigned long int this_block = 0;
    unsigned long int result_temp = 0;
    unsigned long int nr_packages = 2;
    unsigned long int my_results[] = {0, 0, 0};

    MPI_Status my_status[2];

    MPI_Recv(my_block,2*SIZE_OF_BLOCK,MPI_UNSIGNED_LONG,0,DATA,MPI_COMM_WORLD,&status);
    MPI_Irecv(&(my_block[2*SIZE_OF_BLOCK]),SIZE_OF_BLOCK,MPI_UNSIGNED_LONG,0,MPI_ANY_TAG,MPI_COMM_WORLD, &(requests[FOR_RECIVE]));
    
    do {
      
      if (nr_packages == 0)
        break;

      result_temp = 0;

      #pragma omp parallel
      {
        unsigned long int temp_result = 0;

        #pragma omp for
        for (unsigned long int i = 0; i < SIZE_OF_BLOCK; i++) {
          temp_result += checkIfPrime(my_block[i + this_block*SIZE_OF_BLOCK]);
        }

        #pragma omp critical
        {
          result_temp += temp_result;
        }
      }
      // result_temp = 0;
      // for (unsigned long int i = 0; i<SIZE_OF_BLOCK; i++) {
      //   if (checkIfPrime(my_block[i + this_block*SIZE_OF_BLOCK])){
      //     // printf("%ld\n", my_block[i + this_block*SIZE_OF_BLOCK]);
      //     result_temp += 1;
      //   }
        
      // }
      my_results[this_block] = result_temp;
      nr_packages -= 1;

      if (nr_packages == 0)
        break;
      my_recive_tag = false;
      MPI_Waitall(2, requests,my_status);
      // printf("my tag %d\n",my_status[FOR_RECIVE].MPI_TAG);
      
      if (my_status[FOR_RECIVE].MPI_TAG == DATA) {

        nr_packages += 1;
        MPI_Irecv(&(my_block[this_block*SIZE_OF_BLOCK]),SIZE_OF_BLOCK,MPI_UNSIGNED_LONG,0,MPI_ANY_TAG,MPI_COMM_WORLD, &(requests[FOR_RECIVE]));
        my_recive_tag = true;
      }

      MPI_Isend(&(my_results[this_block]),1,MPI_UNSIGNED_LONG,0,RESULT,MPI_COMM_WORLD, &(requests[FOR_SEND]));
      
      this_block = next_block(this_block);
    } while(true);
    MPI_Wait (&(requests[FOR_SEND]), MPI_STATUS_IGNORE);
    MPI_Isend(my_results + this_block,1,MPI_UNSIGNED_LONG,0,RESULT,MPI_COMM_WORLD, &(requests[FOR_SEND]));

    MPI_Wait (&(requests[FOR_SEND]), MPI_STATUS_IGNORE);

    
    MPI_Send(my_results,1,MPI_UNSIGNED_LONG,0,ENDING,MPI_COMM_WORLD);
    printf("slave: ABCD");
  }
  // synchronize/finalize your computations
  if (myrank == 0) {
    printf("Liczb pierwszych: %ld\n", my_result);
    free(numbers);
    free(my_flags_to_not_send);
    free(result_temp);
  }
  else {
    free(my_block);
    free(requests);
  }

  if (!myrank) {
    gettimeofday(&ins__tstop, NULL);
    ins__printtime(&ins__tstart, &ins__tstop, ins__args.marker);
  }
    
  MPI_Finalize();
  
}
