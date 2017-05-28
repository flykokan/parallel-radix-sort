#include <stdio.h>  /* for standard input-output */
#include <stdlib.h> /* for malloc,free,atoi etc... */
#include <mpi.h> /* for MPI standard functions */
#include <omp.h> /* for OMP functions */
#include "MyMPI_OMP.h" /* for my functions */

int main(int argc,char **argv)
{
    int N;   /* the number of keys for sorting */
    int nbuckets; /* # of buckects, one for each possible value of the digit */ 
    int ndigits;   /* # of digits each key will be of, user defines it */ 
    int digit_size;  /* size of each digit in bits */
    int *LCacc,*S,*D;  /* Vectors: S has local keys,LCacc is the local accumulation vector, 
					D is used for the movement step */
	double elapsed_time;  /* counts time */
	int row_start,col_start,split_start; /* indexes of the bounds of process part in global counter Gl */
	int row_end,col_end,split_end;       /* related to the receiving information for that process */
	int *Gl;           /* global counter matrix */
	int *Gl_portion;         /* portion of Gl for this process */
	int **local_Gl_portion;  /* global counter matrix for this process */
    int id;   /* Process ID number */
	int p;    /* Number of processes */ 
	int nthreads;  /* Number of threads */
	int tid;    /* Thread id */
	int i,value,dig,w;
	int *count,*disp,*recv_count,*recv_disp,*send_count,*send_disp; /* count and
																	displacement arrays */
	int local_keys;    /* # of keys held by this process */ 

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&id);
	MPI_Comm_size(MPI_COMM_WORLD,&p);

	if(argc!=4)
	{
		if(!id)
			printf("Usage: ./fname_exe <number_of_digits> <number_of_threads> <input_filename>\n");
		MPI_Finalize();
		exit(1);
	}

	read_keys(argv[3],id,p,&N,&local_keys,&S);
	/* init */
	ndigits=atoi(argv[1]);
	nthreads=atoi(argv[2]);
	omp_set_num_threads(nthreads); 
    digit_size=(sizeof(int)*8)/ndigits; 
    nbuckets=power(2,digit_size);   
	Gl=my_malloc(id,p*nbuckets*(sizeof(int)));
	Gl_portion=my_malloc(id,nbuckets*(sizeof(int)));
	D=my_malloc(id,local_keys*sizeof(int));
	local_Gl_portion=(int **)my_malloc(id,nbuckets*(sizeof(int*))); 
    for(i=0;i<nbuckets;i++)
		local_Gl_portion[i]=my_malloc(id,nthreads*(sizeof(int)));
    /* Start Timer */
    MPI_Barrier(MPI_COMM_WORLD);
	elapsed_time-=MPI_Wtime();
	/******* HERE starts the MPI&OpenMP Radix sort *******/
	#pragma omp parallel private(LCacc,tid,dig,w)
	{
		LCacc=my_malloc(id,nbuckets*sizeof(int));
		tid=omp_get_thread_num();
        for(dig=0;dig<ndigits;dig++)
		{
		/* for least significant digit to most significant digit */
			/* Initialization of counters in the global counter matrix Gl */
            for(w=0;w<nbuckets;w++)
                local_Gl_portion[w][tid]=0;
			/* Compute a histogram of the number of keys for each
             * possible value of the current significant digit */
            #pragma omp for private (i, value) schedule (static)
            for(i=0;i<local_keys;i++)
            {
                value=get_digit_value(S[i],digit_size,dig);
                local_Gl_portion[value][tid]++;
            }
	
            partial_sum(local_Gl_portion,LCacc,nbuckets,tid,nthreads);
			/* Move keys from its part of S to D using its local counters LCacc */
            #pragma omp for private (i, value) schedule (static)
            for(i=0;i<local_keys;i++)
		    {
                value=get_digit_value(S[i],digit_size,dig);
                D[LCacc[value]]=S[i];
                LCacc[value]++;
            }
			/* Compute Gl_portion of this process gathered from threads */
			#pragma omp for private (i,w)
			for(i=0;i<nbuckets;i++)
			{
				Gl_portion[i]=0;
				for(w=0;w<nthreads;w++)
					Gl_portion[i]+=local_Gl_portion[i][w];
			}
			/* there is no data nor functional parallelism */
			/* to exploit for the following functions */ 
            #pragma omp single
			{
			    MPI_Allgather(Gl_portion,nbuckets,MPI_INT,Gl,nbuckets,MPI_INT,MPI_COMM_WORLD);
	
		        counts_displs(id,p,N,&send_count,&send_disp,&recv_count,&recv_disp,Gl,nbuckets,
					&row_start,&col_start,&split_start,&row_end,&col_end,&split_end);
		        MPI_Alltoallv(D,send_count,send_disp,MPI_INT,S,recv_count,recv_disp,MPI_INT,MPI_COMM_WORLD);

		        redistribute_keys(id,p,N,Gl,nbuckets,S,recv_count,row_start,
				    col_start,split_start,row_end,col_end,split_end);
			}
		}
	}
	/******* HERE stops the MPI&OpenMP Radix sort *******/
	/* Stop Timer */
	elapsed_time+=MPI_Wtime();

	print_sorted_keys(id,p,N,S,local_keys);

	if(!id)
		printf("Total elapsed time: %10.6f\n",elapsed_time);

	MPI_Finalize();
    return 0;
}








