#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

int ProcNum;
int ProcRank;


int main(int argc, char* argv[]) {
    double* pMatrix;
    double* pVector;
    double* pResult;
    int Size;
    double* pProcRows;
    double* pProcResult;
    int RowNum;
    double Start, Finish, Duration;
    int i, j;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
//-----------------------------------------------------------------------------------------
    int RestRows;
    if (ProcRank == 0) {
        do {
            Size = atoi(argv[1]);
            if (Size < ProcNum) {
                printf("Size of the objects must be greater than number of processes! \n ");
            }

        }
        while (Size < ProcNum);
    }
    MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    RestRows = Size;
    for (i=0; i<ProcRank; i++)  RestRows = RestRows-RestRows/(ProcNum-i);
    RowNum = RestRows/(ProcNum-ProcRank);
    pVector = (double *) malloc(Size * sizeof(double));
    pResult = (double *) malloc(Size * sizeof(double));
    pProcRows = (double *) malloc(RowNum * Size * sizeof(double));
    pProcResult = (double *) malloc(RowNum * sizeof(double));
    if (ProcRank == 0) {
        pMatrix = (double *) malloc(Size * Size * sizeof(double));
        for (i=0; i<Size; i++) {
            pVector[i] = rand();
            for (j=0; j<Size; j++)
                pMatrix[i*Size+j] = rand();
        }
    }

    Start = MPI_Wtime();
//-----------------------------------------------------------------------
    int *pSendNum;
    int *pSendInd; 
    RestRows=Size;
    MPI_Bcast(pVector, Size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    pSendInd = (int *) malloc(ProcNum * sizeof(int));
    pSendNum = (int *) malloc(ProcNum * sizeof(int));
    RowNum = (Size/ProcNum);
    pSendNum[0] = RowNum*Size;
    pSendInd[0] = 0;
    for (i=1; i<ProcNum; i++) {
      RestRows -= RowNum;
      RowNum = RestRows/(ProcNum-i);
      pSendNum[i] = RowNum*Size;
      pSendInd[i] = pSendInd[i-1]+pSendNum[i-1];
    }
    MPI_Scatterv(pMatrix, pSendNum, pSendInd, MPI_DOUBLE, pProcRows,
      pSendNum[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    free(pSendNum);
    free(pSendInd);
//---------------------------------------------------------------------
    for (i=0; i<RowNum; i++) {
        pProcResult[i] = 0;
        for (j=0; j<Size; j++) pProcResult[i] += pProcRows[i*Size+j]*pVector[j];
    }
//---------------------------------------------------------------------

    int *pReceiveNum;
    int *pReceiveInd; 
    RestRows=Size; 
    pReceiveNum = (int *) malloc(ProcNum * sizeof(int));
    pReceiveInd = (int *) malloc(ProcNum * sizeof(int));

    pReceiveInd[0] = 0;
    pReceiveNum[0] = Size/ProcNum;
    for (i=1; i<ProcNum; i++) {
        RestRows -= pReceiveNum[i-1];
        pReceiveNum[i] = RestRows/(ProcNum-i);
        pReceiveInd[i] = pReceiveInd[i-1]+pReceiveNum[i-1];
    }
    MPI_Allgatherv(pProcResult, pReceiveNum[ProcRank], MPI_DOUBLE, pResult, pReceiveNum, pReceiveInd, MPI_DOUBLE, MPI_COMM_WORLD);
    free(pReceiveNum);
    free(pReceiveInd);
//--------------------------------------------------------------------
    Finish = MPI_Wtime();
    Duration = Finish - Start;
    if (ProcRank == 0) {
        printf("Time of execution = %f\n", Duration);
    }
//--------------------------------------------------------------------
    if (ProcRank == 0) free(pMatrix);
    free(b);
    free(pResult);
    free(pProcRows);
    free(pProcResult);
//--------------------------------------------------------------------
    MPI_Finalize();
    return 0;
}
