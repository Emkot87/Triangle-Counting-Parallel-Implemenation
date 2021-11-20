// Kotoulas Emmanouil 9697 emmakoto@ece.auth.gr Parallel Systems Triangle Counting in adjecency matrices
// compile with gcc pThrdsForBench.c mmio.c -O4 -pthread -o pThrdsForBenchAll.out
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "mmio.h"
#include <pthread.h>
#include <sys/time.h> 
#include <math.h>

// variables to calculate time
struct timeval startwtime, endwtime;
double seq_time;

//defining mutex
pthread_mutex_t dp_mtx;
double dp;
uint32_t counter = 0;

// Structs to contain the various forms of matrices
typedef struct Csr{
    uint32_t n_rows; 
    uint32_t n_cols;
    uint32_t n_nz;
    uint32_t* row_ptrs;
    uint32_t* col_indices;
    uint32_t* final_vals;
    uint32_t  iter;
    uint32_t  threads;
}Csr;

typedef struct Coo{
    uint32_t n_rows;
    uint32_t n_cols;
    uint32_t n_nz;
    uint32_t* row_indices;
    uint32_t* col_indices;
}Coo;

// Function reading mtx file in coo format and passing it to a coo struct
Coo MtxToCoo(int argc,char *argv){
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    uint32_t M, N, nz;   
    uint32_t i, *I, *J;

    if ((f = fopen(argv, "r")) == NULL) 
        exit(1);
    

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }


    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        exit(1);

    /* reseve memory for matrices */

    I = (uint32_t *) malloc(nz * sizeof(uint32_t));
    J = (uint32_t *) malloc(nz * sizeof(uint32_t));


    for (i=0; i<nz; i++)
    {
        fscanf(f, "%d %d \n", &I[i], &J[i] );
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
    }

    if (f !=stdin) fclose(f);

    Coo sparseMatrix ;
    sparseMatrix.n_rows = M ;
    sparseMatrix.n_cols = M ;
    sparseMatrix.n_nz = nz ;
    sparseMatrix.col_indices = J;
    sparseMatrix.row_indices = I;   

    return sparseMatrix;
}

// Function by Dimitris Floros to make a csc from a coo
void coo2csc(
  uint32_t       * const row,       /*!< CSC row start indices */
  uint32_t       * const col,       /*!< CSC column indices */
  uint32_t const * const row_coo,   /*!< COO row indices */
  uint32_t const * const col_coo,   /*!< COO column indices */
  uint32_t const         nnz,       /*!< Number of nonzero elements */
  uint32_t const         n,         /*!< Number of rows/columns */
  uint32_t const         isOneBased /*!< Whether COO is 0- or 1-based */
) {

  // ----- cannot assume that input is already 0!
  for (uint32_t l = 0; l < n+1; l++) col[l] = 0;


  // ----- find the correct column sizes
  for (uint32_t l = 0; l < nnz; l++)
    col[col_coo[l] - isOneBased]++;

  // ----- cumulative sum
  for (uint32_t i = 0, cumsum = 0; i < n; i++) {
    uint32_t temp = col[i];
    col[i] = cumsum;
    cumsum += temp;
  }
  col[n] = nnz;
  // ----- copy the row indices to the correct place
  for (uint32_t l = 0; l < nnz; l++) {
    uint32_t col_l;
    col_l = col_coo[l] - isOneBased;

    uint32_t dst = col[col_l];
    row[dst] = row_coo[l] - isOneBased;

    col[col_l]++;
  }
  // ----- revert the column pointers
  for (uint32_t i = 0, last = 0; i < n; i++) {
    uint32_t temp = col[i];
    col[i] = last;
    last = temp;
  }

}

// Triangle calculating function takes a pointer to struct and calculates the product given the number of threads
void* TriangleCounting(void * arg){
    Csr *csrSparse = (Csr*)arg;
    uint32_t start , end ;
    start = (csrSparse->iter)*((csrSparse->n_cols + 1)/csrSparse->threads) ;
    end = (csrSparse->iter + 1)*((csrSparse->n_cols + 1)/csrSparse->threads);

    
    if(csrSparse->iter == (csrSparse->threads - 1)){
        end = csrSparse->n_cols + 1 ;
    }


    for( uint32_t g = start ; g < end ; g++ ){
        for( uint32_t j = csrSparse->row_ptrs[g] ; j < csrSparse->row_ptrs[g+1] ; j++ ){
            for(uint32_t k = csrSparse->row_ptrs[csrSparse->col_indices[j]] ; k < csrSparse->row_ptrs[csrSparse->col_indices[j]+1] ; k++){
                for(uint32_t l = csrSparse->row_ptrs[g] ; l < csrSparse->row_ptrs[g+1] ; l++ ){
                    if(csrSparse->col_indices[k] == csrSparse->col_indices[l]){
                        pthread_mutex_lock(&dp_mtx);

                            csrSparse->final_vals[j]++;

                        pthread_mutex_unlock(&dp_mtx);
                    }
                }
            }
        }
    }
    pthread_exit(NULL);
}


int main(int argc, char *argv[]){ 
    uint32_t i,j,k,g,l,test;

    // read mtx
    char matrices[5][20] ={"belgium_osm.mtx","com-Youtube.mtx","dblp-2010.mtx","mycielskian13.mtx","NACA0015.mtx"};
    printf("time,sum,num of triangles,threads,matrix\n");
    for(test = 0 ; test <5 ; test++){

        Coo cooHalf = MtxToCoo(argc,matrices[test]);

        //allocate space for the whole matrix
        uint32_t * cooWhole_row = (uint32_t *)malloc(cooHalf.n_nz * 2 * sizeof(uint32_t));
        uint32_t * cooWhole_col = (uint32_t *)malloc(cooHalf.n_nz * 2 * sizeof(uint32_t));


        // pass the indices and add them to the other array to have the whole matrix in an unordered way

        for( i = 0 ; i < cooHalf.n_nz ; i++){
            cooWhole_row[i] = cooHalf.row_indices[i];
            cooWhole_row[i + cooHalf.n_nz] = cooHalf.col_indices[i];
            cooWhole_col[i] = cooHalf.col_indices[i];
            cooWhole_col[i + cooHalf.n_nz] = cooHalf.row_indices[i];
        }

        //make a new struct and pass the values and the indices  
        Coo cooWhole ;
        cooWhole.row_indices = cooWhole_row ;
        cooWhole.col_indices = cooWhole_col ;
        cooWhole.n_cols = cooHalf.n_cols;
        cooWhole.n_nz = 2*cooHalf.n_nz ;

        // allocate space for the csc arrays
        uint32_t * csc_row = (uint32_t *)malloc(cooWhole.n_nz  * sizeof(uint32_t));
        uint32_t * csc_col = (uint32_t *)malloc((cooWhole.n_cols + 1) * sizeof(uint32_t));

        // fill the csc arrays from the coo
        coo2csc(csc_row, csc_col,cooWhole.col_indices,cooWhole.row_indices,cooWhole.n_nz,cooWhole.n_cols,0);
        
        // table containing various thread values to benchmark with
        int threadsTable[9] = {1,2,4,8,16,32,64,256,1024};
        for(int thrd = 0 ; thrd < 9 ; thrd++){
            for(int many = 0 ; many < 3 ; many++ ){
                int n = threadsTable[thrd];

                uint32_t * productValues = (uint32_t *)malloc(cooWhole.n_nz * sizeof(uint32_t));

                // make an array csr (which is equal to a csc in case of symmetric matrices) structs and pass the csc arrays
                Csr *csrSparse = malloc (n*sizeof(Csr));

                for(int i = 0 ; i < n ; i++){
                    csrSparse[i].row_ptrs = csc_col ;
                    csrSparse[i].col_indices = csc_row ;
                    csrSparse[i].iter = 0 ;
                    csrSparse[i].n_cols = cooWhole.n_cols;
                    csrSparse[i].final_vals = productValues;
                    csrSparse[i].threads = threadsTable[thrd];
                }

                //Initialize mutex and thread array
                pthread_attr_t pthread_custom_attr;
                dp = 0.0;
                pthread_mutex_init (&dp_mtx, NULL);
                pthread_t *threads;
                threads = (pthread_t*)malloc(n*sizeof(pthread_t));

                
                for (uint32_t l = 0; l < cooWhole.n_nz ; l++) 
                    productValues[l] = 0 ;
                //get starting time
                gettimeofday (&startwtime, NULL);
                
                //create the threads and call the function, pass i to the current csr struct
                for (i = 0; i < n; i++) {
                    csrSparse[i].iter = i;
                    pthread_create(&threads[i] ,NULL ,TriangleCounting ,(void *)(csrSparse+i));
                }
                
                //call to join the threads
                for (i = 0; i<n; i++) {
                    pthread_join(threads[i],NULL);
                }

                //free the thread matrix
                free(threads);

                pthread_mutex_destroy(&dp_mtx);
            
                uint32_t cumsum = 0;
                // sum the values of the product of the multiplication
                for(j = 0 ; j < cooWhole.n_nz ; j++){
                    cumsum += productValues[j];
                }

                // calculate triangles
                uint32_t triangles = cumsum/6;
                // get final time
                gettimeofday (&endwtime, NULL);
                seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
                //print final results in csv format
                printf("%f,%u,%u,%d,%s\n",seq_time,cumsum, triangles, threadsTable[thrd],matrices[test]);
            }
        }
    }
}
