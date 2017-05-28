#include <stdio.h>  /* for standard input-output */
#include <stdlib.h> /* for malloc,free,atoi etc... */
#include <time.h> /* for time,gettimeofday */
#include <sys/time.h> /* for struct timeval */ 

int get_digit_value(int num,int digit_size,int dig);
void swap_addresses(int **S,int **D);
int power(int base,int e);

int main(int argc,char **argv)
{
    int N;   /* the number of keys for sorting */
    int nbuckets; /* # of buckects, one for each possible value of the digit */ 
    int ndigits;   /* # of digits each key will be of, user defines it */ 
    int digit_size;  /* size of each digit in bits */
    int *S,*C,*D; /* Vectors: S has the keys, C has the istogram with possible
                values of the current digit, D is used for the movement step */
    int tmp,accum,i,value,dig;
    FILE *fin,*fout;   /* pointers to input and output file */
    struct timeval start,stop;  /* for gettimeofday */

    if(argc!=3)
    {
        printf("Usage: ./fname.exe <number_of_digits> <input_filename>\n");
        exit(1);
    }
       
	if((fin=fopen(argv[2],"r"))==NULL)
    {
        printf("Cannot open file %s .\n",argv[2]);
        exit(1);
	}
    if((fout=fopen("seq_radix.txt","w"))==NULL)
    {
        printf("Cannot open file seq_radix.txt.\n");
        exit(1);
	}

	fscanf(fin,"%d",&N); /* read the number of keys */
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
    digit_size=(sizeof(int)*8)/ndigits;   /* compute digit_size, 
									         ndigits is known */
    nbuckets=power(2,digit_size);   /* 2^digit_size buckets needed */
    D=(int *)malloc(N*sizeof(int));

    gettimeofday(&start,NULL);    
/******* HERE starts the Sequential Radix sort *******/
    for(dig=0;dig<ndigits;dig++)
    {
    /* for least significant digit to most significant digit */
        C=(int *)malloc(nbuckets*sizeof(int));  /* initialize C */
		if(!C)
		{
			printf("Allocation failed");
		    exit(1);
	    } 
        /*** This is the count step ***/ 
        for(i=0;i<N;i++)
        {
            value=get_digit_value(S[i],digit_size,dig);
            C[value]++;
        }
        /*** This is the accummulation step ***/
        tmp=C[0];
        C[0]=0;
        for(i=1;i<nbuckets;i++)
        {
            accum=tmp+C[i-1];
            tmp=C[i];
            C[i]=accum;
        }
        /*** This is the movement step ***/
        for(i=0;i<N;i++)
        {
            value=get_digit_value(S[i],digit_size,dig);
            D[C[value]]=S[i];
            C[value]++;
       }
       swap_addresses(&S,&D);  /* S,D exchange their role */
    }
/******* HERE stops the Sequential Radix sort *******/
    gettimeofday(&stop,NULL);
    
    printf("Time of Sequential Radix sort is %ld s",(stop.tv_sec-start.tv_sec));
    printf("  %ld micros\n",(abs(stop.tv_usec-start.tv_usec)));
	/* write sorted keys to seq_radix.txt */
    for(i=0;i<N;i++)
        fprintf(fout,"%d ",S[i]);
    free(D);
    free(C);
    free(S);
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


