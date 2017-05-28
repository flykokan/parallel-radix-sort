/*
 * create_input.c -- Using time for seed of rand 
 * a file is created consisting of an integer,
 * the size of the problem,followed by that number 
 * of random integers
 */

#include <stdio.h>  /* for standard input-output */
#include <stdlib.h> /* for malloc,free,atoi etc... */

int main(int argc,char **argv)
{
	int utime; 
    long ltime;
	int N,*S,i;
	FILE *fp;

	if(argc!=3)
    {
        printf("Usage: ./fname.exe <size_problem(N)> <output_filename>\n");
        exit(1);
    }
	N=atoi(argv[1]);
    ltime=time(NULL);  
    utime=(unsigned int) ltime;
    srand(utime);
    S=(int *)malloc(N*sizeof(int));
	if(!S)
	{
		printf("Allocation failed");
		exit(1);
	}
	if((fp=fopen(argv[2],"w"))==NULL)
    {
        printf("Cannot open file.\n");
        exit(1);
    }
    fprintf(fp,"%d ",N);
    for(i=0;i<N;i++)
	{
        *(S+i)=rand();
		fprintf(fp,"%d ",*(S+i));
	}
	free(S);
	fclose(fp);
	return 0;
}
