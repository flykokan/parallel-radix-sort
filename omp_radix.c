#include <stdio.h>  /* for standard input-output */
#include <stdlib.h> /* for malloc,free,atoi etc... */
#include <time.h> /* for time,gettimeofday */
#include <sys/time.h> /* for struct timeval */
#include <omp.h>  /* for OpenMP functions */ 

int get_digit_value(int num,int digit_size,int dig);
void partial_sum(int **Gl,int *LCacc,int nbuckets,int tid,int nthreads);
void swap_addresses(int **S,int **D);
int power(int base,int e); 

int main(int argc,char **argv)
{
    int N;   /* the number of keys for sorting */
    int nbuckets; /* # of buckects, one for each possible value of the digit */ 
    int ndigits;   /* # of digits each key will be of, user defines it */ 
    int digit_size;  /* size of each digit in bits */
    int *S,*LCacc,*D; /* Vectors: S has the keys,LCacc is the local accumulation vector, 
					  D is used for the movement step */
	int **Gl;   /* Global counter matrix */
    int tid;   /* Thread id */
	int nthreads;  /* Number of threads */  
	int i,w,value,dig;
    FILE *fin,*fout;   /* pointers to input and output file */
    struct timeval start,stop;  /* for gettimeofday */

    if(argc!=4)
    {
        printf("Usage: ./fname.exe <number_of_digits> <number_of_threads> <input_filename>\n");
        exit(1);
    }
       
	if((fin=fopen(argv[3],"r"))==NULL)
    {
        printf("Cannot open file %s .\n",argv[3]);
        exit(1);
	}
    if((fout=fopen("omp_radix.txt","w"))==NULL)
    {
        printf("Cannot open file omp_radix.txt.\n");
        exit(1);
	}
    
	fscanf(fin,"%d",&N);  /* read the number of keys */
    S=(int *)malloc(N*sizeof(int));
	D=(int *)malloc(N*sizeof(int));
	if(!S || !D)
	{
		printf("Allocation failed");
		exit(1);
	}
	/* read the keys */
	for(i=0;i<N;i++)
        fscanf(fin,"%d ",S+i);

    ndigits=atoi(argv[1]);     
	nthreads=atoi(argv[2]);
	omp_set_num_threads(nthreads);  /* define the # of threads */
    digit_size=(sizeof(int)*8)/ndigits;   /* compute digit_size, 
									         ndigits is known */
    nbuckets=power(2,digit_size);   /* 2^digit_size buckets needed */
	/* init */
	LCacc=(int *)malloc(nbuckets*sizeof(int));
	if(!LCacc)
	{
		printf("Allocation failed");
		exit(1);
	}
   	Gl=(int **)malloc(nbuckets*(sizeof(int*))); 
    for(i=0;i<nbuckets;i++)
		Gl[i]=(int *)malloc(nthreads*(sizeof(int)));

    gettimeofday(&start,NULL);
	/******* HERE starts the Parallel OpenMP Radix sort *******/
    #pragma omp parallel private(LCacc,tid,dig,w)
	{
		LCacc=(int *)malloc(nbuckets*sizeof(int));
		tid=omp_get_thread_num();
        for(dig=0;dig<ndigits;dig++)
		{
		/* for least significant digit to most significant digit */
			/* Initialization of counters in the global counter matrix Gl */
            for(w=0;w<nbuckets;w++)
                Gl[w][tid]=0;
			/* Compute a histogram of the number of keys for each
             * possible value of the current significant digit */
            #pragma omp for private (i, value) schedule (static)
            for(i=0;i<N;i++)
            {
                value=get_digit_value(S[i],digit_size,dig);
                Gl[value][tid]++;
            }
	
            partial_sum(Gl,LCacc,nbuckets,tid,nthreads);
			/* Move keys from its part of S to D using its local counters LCacc */
            #pragma omp for private (i, value) schedule (static)
            for(i=0;i<N;i++)
		    {
                value=get_digit_value(S[i],digit_size,dig);
                D[LCacc[value]]=S[i];
                LCacc[value]++;
            }
			/* One of the processors exchanges the roles of vectors S and D */
            #pragma omp single
            swap_addresses(&S,&D);
	    }
	}
/******* HERE stops the Parallel OpenMP Radix sort *******/
    gettimeofday(&stop,NULL);
    
    printf("Time of OMP Radix sort is %ld s",(stop.tv_sec-start.tv_sec));
    printf("  %ld micros\n",(abs(stop.tv_usec-start.tv_usec)));
	/* write sorted keys to omp_radix.txt */
	for(i=0;i<N;i++)
        fprintf(fout,"%d ",S[i]);

    free(D);
    free(S);
    free(Gl);
	free(LCacc);
    fclose(fin);
	fclose(fout);
    return 0;
}

/* Given  the size of digit(digit_size) and the
digit No(dig), this function returns the value(for
vector C) of the integer num of current dig */
   
int get_digit_value(int num,int digit_size,int dig)
{
    int value,mask;
    
    value=num >> (digit_size*dig);
    mask=power(2,digit_size)-1;
    value=value & mask;
    return value;
}

void swap_addresses(int **S,int **D)
{
    int *temp;
    
    temp=*S;
    *S=*D;
    *D=temp;
}


int power(int base,int e)
{
    int temp;
    
    temp=1;
    for(;e>0;e--)
       temp=temp*base;
       
    return temp;
}

/* This function computes the local accumulation vector LCacc using 
 * global counters Gl.Processor Ptid knows the first place where the first
 * key of each of its buckets should be written. So, the first key belonging
 * to bucket i of processor Ptid is placed in the position of vector D just
 * after all the keys belonging to buckets 0 to i-1, and after all the keys
 * of the same bucket i belonging to processors P0 to Ptid-1
 */

void partial_sum(int **Gl,int *LCacc,int nbuckets,int tid,int nthreads)
{
	int sum1,sum2;
	int k,l,ip;

	for(k=0;k<nbuckets;k++)
	{
		sum1=0;
		sum2=0;
	    for(ip=0;ip<nthreads;ip++)
	    {
			for(l=0;l<k;l++)
				sum1+=Gl[l][ip];
		}
		for(ip=0;ip<tid;ip++)
			sum2+=Gl[k][ip];
		LCacc[k]=sum1+sum2;
	}
}
