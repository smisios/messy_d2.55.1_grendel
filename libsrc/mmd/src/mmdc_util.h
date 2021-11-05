#ifndef MMD_FEFS_H
#define MMD_FEFS_H 1

#define my_int  int
#define my_PINT PINT
#define my_char STRING

#define MMD_MAX_MODEL 64

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <mpi.h>

#include "cfortran.h"
struct BufDef {
  int                  BufLen;       /* Length of buffer                    */
  int                  pBufLen;      /* Length of parent buffer             */
  long                 BufIndex;     /* index in Send Buffer                */ 
  double               *SendBuf;     /* Pointer of Data in Send buffer      */
 
  struct BufDef      *next; 
};

struct ModelDef {
  MPI_Aint             TotalBufferSize;
  MPI_Comm             model_comm;    /* Communicator of this model          */
  MPI_Comm             inter_comm;    /* Inter communicator model and client */
  MPI_Comm             intra_comm;    /* Intra communicator model and client */
  int                  model_rank;    /* Rank of this model                  */
  int                  model_npes;    /* Number of PEs this model            */
  int                  inter_npes;    /* Number of PEs of remote model       */
  MPI_Win              BufWin;        /* MPI RMA windows                     */
  struct BufDef        *buf;
};

double MMDc_U_Time (void);

#endif
