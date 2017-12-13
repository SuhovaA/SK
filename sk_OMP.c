#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#ifdef _OPENMP
    #include <omp.h>
#endif
void prt1a(char *t1, double *v, int n,char *t2) ;
void wtime(double *t)
{
  static int sec = -1;
  struct timeval tv;
  gettimeofday(&tv, 0);
  if (sec < 0) sec = tv.tv_sec;
  *t = (tv.tv_sec - sec) + 1.0e-6*tv.tv_usec;
}

int N;
int i, j, k;
double *A;
#define A(i,j) A[(i)*(N)+(j)]
double *X, *b;
int NUM_THREADS = 2;
int main(int argc,char **argv){
        double time0, time1;

        if (argc > 1)
        {
            NUM_THREADS = atoi(argv[1]);
            N = atoi(argv[2]);
        }
        /*FILE *in;
        in=fopen("data.in","r");
        if(in==NULL) {
                printf("Can not open 'data.in' "); exit(1);
        }

        i=fscanf(in,"%d", &N);
        if(i<1) {
                printf("Wrong 'data.in' (N ...)"); exit(2);
        }
        */
        A = (double *) malloc(N * N * sizeof(double));
        b = (double *) malloc(N * sizeof(double));
        X = (double *) malloc(N * sizeof(double));
        printf("A x b ( N = %d )\n----------------------------------\n",N);
        wtime(&time0);
        /* initialize array A*/
        for (i = 0; i <= N - 1; i++) {
            for(j = 0; j <= N - 1; j++) {
                A(i,j) = rand();
                //printf("%lf ", A[i * N + j]);
            }
            //printf("\n");
        }
        for (i = 0; i < N; i++) {
            b[i] = rand();
            //printf("%lf ", b[i]);
        }
        printf("\n");
        int sum;
#ifdef _OPENMP
        omp_set_num_threads(NUM_THREADS);
#endif
        #pragma omp parallel for private(j, sum) schedule(dynamic) num_threads(NUM_THREADS)
        for(i = 0; i < N; i++) {
                sum = 0;
                for(j = 0; j < N; j++) {
                        sum += A(i,j) * b[j];
                }
                X[i] = sum;
        }

        wtime(&time1);
#ifdef _OPENMP
        printf("Number of threads=%d %d\n",omp_get_max_threads(),_OPENMP);
#endif
        printf("Time in seconds=%gs\n",time1-time0);
        prt1a("X = ( ", X,N>9?9:N,")\n");
        free(A);
        free(b);
        free(X);
        //free(X);
        return 0;
}

void prt1a(char * t1, double *v, int n,char *t2){
        int j;
        printf("%s",t1);
        for(j=0;j<n;j++)
            printf("%.4g%s",v[j], j%10==9? "\n": " ");
        printf("%s",t2);
}
