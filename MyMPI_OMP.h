/**************MACROS******************/
#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) (BLOCK_HIGH(id,p,n)- BLOCK_LOW(id,p,n)+1)

/*******************MISCELLANEOUS FUNCTIONS***************/
void redistribute_keys(int id,int p,int n,int *Gl,int nbuckets,int *S,int *recv_counts,int row_start,
					   int col_start,int split_start,int row_end,int col_end,int split_end);
int get_digit_value(int num,int digit_size,int dig);
int power(int base,int e);
int *my_malloc(int id,int bytes);
void partial_sum(int **Gl,int *LCacc,int nbuckets,int tid,int nthreads);

/**************DATA DISTRIBUTION FUNCTIONS****************/
void counts_displacements(int id,int p,int n,int **count,int **disp);
void counts_displs(int id,int p,int n,int **scounts,int **sdispls,int **rcounts,int **rdispls,int *Gl,int nbuckets,
				   int *row_start,int *col_start,int *split_start,int *row_end,int *col_end,int *split_end);

/*****************INPUT FUNCTIONS**************/
void read_keys(char *in_fname,int id,int p,int *N,int *local_keys,int **S);

/*****************OUTPUT FUNCTIONS*************/
void print_sorted_keys(int id,int p,int N,int *S,int local_keys);
