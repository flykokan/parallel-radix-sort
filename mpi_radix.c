#include <stdio.h>  /* for standard input-output */
#include <stdlib.h> /* for malloc,free,atoi etc... */
#include <mpi.h> /* for MPI standard functions */
#include "MyMPI.h" /* for my functions */

int main(int argc,char **argv)
{
    int N;   /* the number of keys for sorting */
    int nbuckets; /* # of buckects, one for each possible value of the digit */ 
    int ndigits;   /* # of digits each key will be of, user defines it */ 
    int digit_size;  /* size of each digit in bits */
    int *LC,*S,*D;  /* Vectors: S has local keys,LC is the local counter vector, 
					D is used for the movement step */
	double elapsed_time;  /* counts time */
	int row_start,col_start,split_start; /* indexes of the bounds of process part in global counter Gl */
	int row_end,col_end,split_end;       /* related to the receiving information for that process */
	int *Gl;           /* global counter matrix */
	int *Gl_portion;         /* portion of Gl for this process */  
    int id;   /* Process ID number */
	int p;    /* Number of processes */ 
	int i,value,dig,tmp,accum;
	int *count,*disp,*recv_count,*recv_disp,*send_count,*send_disp; /* count and
																	displacement arrays */
	int local_keys;    /* # of keys held by this process */ 

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&id);
	MPI_Comm_size(MPI_COMM_WORLD,&p);

	if(argc!=3)
	{
		if(!id)
			printf("Usage: ./fname_exe <number_of_digits> <input_filename>\n");
		MPI_Finalize();
		exit(1);
	}

	read_keys(argv[2],id,p,&N,&local_keys,&S);
    /* init */
	ndigits=atoi(argv[1]);
    digit_size=(sizeof(int)*8)/ndigits; 
    nbuckets=power(2,digit_size);  
	Gl_portion=my_malloc(id,nbuckets*(sizeof(int))); 
	Gl=my_malloc(id,p*nbuckets*(sizeof(int)));
	D=my_malloc(id,local_keys*sizeof(int));

    /* Start Timer */
    MPI_Barrier(MPI_COMM_WORLD);
	elapsed_time-=MPI_Wtime();
	/******* HERE starts the MPI Radix sort *******/
    for(dig=0;dig<ndigits;dig++)
	{
		/* for least significant digit to most significant digit */
        LC=my_malloc(id,nbuckets*sizeof(int));
        for(i=0;i<nbuckets;i++)
                LC[i]=0;
		/*** This is the count step ***/
        for(i=0;i<local_keys;i++)
        {
            value=get_digit_value(S[i],digit_size,dig);
            LC[value]++;
        }
		/*** This is the accummulation step ***/
		tmp=LC[0];
        LC[0]=0;
        for(i=1;i<nbuckets;i++)
        {
            accum=tmp+LC[i-1];
            tmp=LC[i];
            LC[i]=accum;
        }
		/*** This is the movement step ***/
        for(i=0;i<local_keys;i++)
		{
            value=get_digit_value(S[i],digit_size,dig);
            D[LC[value]]=S[i];
            LC[value]++;
		}
		/* Compute local Gl_portion */ 
		Gl_portion[0]=LC[0];
		for(i=1;i<nbuckets;i++)
			Gl_portion[i]=LC[i]-LC[i-1];

	    MPI_Allgather(Gl_portion,nbuckets,MPI_INT,Gl,nbuckets,MPI_INT,MPI_COMM_WORLD);
		
		counts_displs(id,p,N,&send_count,&send_disp,&recv_count,&recv_disp,Gl,nbuckets,
			  &row_start,&col_start,&split_start,&row_end,&col_end,&split_end);
		MPI_Alltoallv(D,send_count,send_disp,MPI_INT,S,recv_count,recv_disp,MPI_INT,MPI_COMM_WORLD);
		
		redistribute_keys(id,p,N,Gl,nbuckets,S,recv_count,row_start,
			col_start,split_start,row_end,col_end,split_end);
	}
	/******* HERE stops the OMP Radix sort *******/
	/* Stop Timer */
	elapsed_time+=MPI_Wtime();

	print_sorted_keys(id,p,N,S,local_keys);

	if(!id)
		printf("Total elapsed time: %10.6f\n",elapsed_time);

	MPI_Finalize();
    return 0;
}
