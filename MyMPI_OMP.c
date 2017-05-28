#include <stdio.h>  /* for standard input-output */
#include <stdlib.h> /* for malloc,free,atoi etc... */
#include <mpi.h> /* for MPI standard functions */
#include "MyMPI_OMP.h" /* for my functions */

/* This function reads the keys from a file.The first element
* of the file is the # of the keys to be sorted followed by the keys 
* themselves.MPI process with id equal to p-1 reads the total keys 
* and distribute them properly to all processes(including itself)
* using MPI function MPI_Scatterv.
*/

void read_keys(char *in_fname,int id,int p,int *N,int *local_keys,int **S)
{
	int *buffer,*count,*disp;
    FILE *fin;
	int i;

    if(id==(p-1))
	{
		if((fin=fopen(in_fname,"r"))==NULL)
			MPI_Abort(MPI_COMM_WORLD,-1);
		else; 
			fscanf(fin,"%d",N);
	}
	MPI_Bcast(N,1,MPI_INT,p-1,MPI_COMM_WORLD);
	*local_keys=BLOCK_SIZE(id,p,*N);
    *S=my_malloc(id,(*local_keys)*sizeof(int));
	if(id==(p-1))
		buffer=my_malloc(id,(*N)*sizeof(int));
	counts_displacements(id,p,*N,&count,&disp);
	if(id==(p-1))
	{
		for(i=0;i<(*N);i++)
			fscanf(fin,"%d",buffer+i);
	}
	MPI_Scatterv(buffer,count,disp,MPI_INT,
		*S,*local_keys,MPI_INT,p-1,MPI_COMM_WORLD);

	free(count);
	free(disp);
	if(id==(p-1))
	{
		free(buffer);
		fclose(fin);
	}
}

/* This function prints in a file the sorted keys after
 * gathering sorted parts S from each process to process
 * with id equal to p-1 using MPI function MPI_Gatherv
 */

void print_sorted_keys(int id,int p,int N,int *S,int local_keys)
{
	int i;
	int *buffer;
	int *count,*disp;
	FILE *fout;

	if(id==(p-1))
		buffer=my_malloc(id,N*sizeof(int));
	counts_displacements(id,p,N,&count,&disp);
	MPI_Gatherv(S,local_keys,MPI_INT,buffer,count,disp,MPI_INT,p-1,MPI_COMM_WORLD);
	if(id==(p-1))
	{
		if((fout=fopen("mpi_omp_radix.txt","w"))==NULL)
			MPI_Abort(MPI_COMM_WORLD,-1);
		else
		{
			for(i=0;i<N;i++)
				fprintf(fout,"%d ",*(buffer+i));
		}
	}
	free(count);
	free(disp);
	if(id==(p-1))
	{
		free(buffer);
		fclose(fout);
	}
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

int power(int base,int e)
{
    int temp;
    
    temp=1;
    for(;e>0;e--)
       temp=temp*base;
       
    return temp;
}

/* This function creates the count and displacement arrays 
 * needed by scatter and gather functions
 */

void counts_displacements(int id,int p,int n,int **count,int **disp)
{
	int i;

	*count=my_malloc(id,p*sizeof(int));
	*disp=my_malloc(id,p*sizeof(int));
	(*count)[0]=BLOCK_SIZE(0,p,n);
	(*disp)[0]=0;
	for(i=1;i<p;i++)
	{
		(*disp)[i]=(*disp)[i-1]+(*count)[i-1];
		(*count)[i]=BLOCK_SIZE(i,p,n);
	}
}

/* This function is called when a process wants to allocate
 * bytes of memory.If the memory allocation fails, the process prints
 * an error message and then aborts execution of the program
 */

int *my_malloc(int id,int bytes)
{
	int *buffer;
	if((buffer=(int *)malloc(bytes*sizeof(int)))==NULL)
	{
		printf("Error:Malloc failed for process %d\n",id);
		fflush(stdout);
		MPI_Abort(MPI_COMM_WORLD,-2);
	}
	return buffer;
}

/* This function creates the count and displacement arrays, 
 * both sending and receiving , needed by our MPI_Alltoallv function.
 * Each process scans the global matrix Gl to get the needed data,
 * starting from bucket of digit 0 of process 0,then the one of digit 0
 * of process 1 etc.It also stores the bounds of process part in 
 * global counter Gl related to the receiving information for that process. 
 * In addition it stores information after splitting buckets where needed in 
 * order to achieve a load balance.
 */

void counts_displs(int id,int p,int n,int **scounts,int **sdispls,int **rcounts,int **rdispls,int *Gl,int nbuckets,
				   int *row_start,int *col_start,int *split_start,int *row_end,int *col_end,int *split_end)
{
	int i,k,j;
	int i_start,j_start;
	int splitting,sum,my_sum;

	*scounts=my_malloc(id,p*sizeof(int));
	*sdispls=my_malloc(id,p*sizeof(int));
	*rcounts=my_malloc(id,p*sizeof(int));
	*rdispls=my_malloc(id,p*sizeof(int));

	my_sum=sum=splitting=0;
	i_start=j_start=i=j=0;
	for(k=0;k<p;k++)
	{
       (*rcounts)[k]=0;
	   (*scounts)[k]=0;
	}
	for(k=0;k<p;k++)
	{
		if((my_sum==(sum-splitting))&&(id==k))
		{
	    /* These actions are for the process part in global counter Gl 
		 * related to the receiving information for that process.  */
			*row_start=j_start;
			*col_start=i_start;
			*split_start=splitting;
			if(splitting>0)
			{
				if(j==0)
					(*rcounts)[p-1]+=splitting;
				else
					(*rcounts)[j-1]+=splitting;
			}
			for(i=i_start;i<nbuckets&&sum<my_sum+BLOCK_SIZE(k,p,n);i++)
			{
				for(j=j_start;j<p&&sum<my_sum+BLOCK_SIZE(k,p,n);j++)
				{
					(*rcounts)[j]+=Gl[j*nbuckets+i];
                    sum+=Gl[j*nbuckets+i];
					if(j==p-1)
						j_start=0;
				}
			}
			splitting=sum-my_sum-BLOCK_SIZE(k,p,n);
			if(splitting>0)
                (*rcounts)[j-1]-=splitting;
			(*scounts)[id]=(*rcounts)[id];
			*row_end=j-1;
			*col_end=i-1;
			*split_end=splitting;
            if(j!=p)
				i--;
			else j=0;
			i_start=i;
			j_start=j;
			
		}
		else
		{
		/* These actions are for the process part in global counter Gl 
		 * related to the sending information for that process.  */
			if(splitting>0)
			{
				if(j==0&&id==p-1)
					
					(*scounts)[k]+=splitting;
				else if(id==j-1)
					(*scounts)[k]+=splitting;
				else;
			}
            for(i=i_start;i<nbuckets&&sum<my_sum+BLOCK_SIZE(k,p,n);i++)
			{
				for(j=j_start;j<p&&sum<my_sum+BLOCK_SIZE(k,p,n);j++)
				{
					if(j==id)
						(*scounts)[k]+=Gl[j*nbuckets+i];
					sum+=Gl[j*nbuckets+i];
					if(j==p-1)
						j_start=0;

				}
			}
			splitting=sum-my_sum-BLOCK_SIZE(k,p,n);
			if(splitting>0)
			{
				if((j-1)==id)
					(*scounts)[k]-=splitting;
			}
			if(j!=p)
				i--;
			else j=0;
			i_start=i;
			j_start=j;
		}
		my_sum+=BLOCK_SIZE(k,p,n);
	}
	(*sdispls)[0]=0;
	(*rdispls)[0]=0;
	for(i=1;i<p;i++)
	{
		(*sdispls)[i]=(*sdispls)[i-1]+(*scounts)[i-1];
		(*rdispls)[i]=(*rdispls)[i-1]+(*rcounts)[i-1];
	}
}

/* This function redistribute the local keys for each process 
 * for the current digit after the alltoall data communication.
 * Using the data from the previous counts_displs function and
 * the global counter Gl, it sorts the received buckets of keys
 */

void redistribute_keys(int id,int p,int n,int *Gl,int nbuckets,int *S,int *recv_counts,int row_start,
					   int col_start,int split_start,int row_end,int col_end,int split_end)
{
	int i,j,k,*ptr;
	int row_split;
	int **chunks;

	ptr=S;
	chunks=(int **)malloc(p*(sizeof(int*)));
	for(i=0;i<p;i++)
		chunks[i]=(int *)malloc(recv_counts[i]*(sizeof(int)));
	for(i=0;i<p;i++)
		for(j=0;j<recv_counts[i];j++)
			*(chunks[i]+j)=*ptr++;
	
	ptr=S;
	Gl[row_end*nbuckets+col_end]-=split_end;
	if(split_start>0)
	{
		if(row_start==0)
			row_split=p-1;
		else 
			row_split=row_start-1;
		for(i=0;i<split_start;i++)
			*ptr++=*chunks[row_split]++;
	}
	for(i=col_start;i<col_end+1;i++)
	{
		for(j=row_start;j<p;j++)
		{
			if(i==col_end&&j>row_end)
				break;
			if(j==p-1)
				row_start=0;
			for(k=0;k<Gl[j*nbuckets+i];k++)
				*ptr++=*chunks[j]++;
		}
	}
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
