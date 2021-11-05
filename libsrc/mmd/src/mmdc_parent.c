 
#include "mmdc_util.h"

/* ------------------------------------------------------------------------*/
/* VARIABLE DECLARATION                                                    */
/* ------------------------------------------------------------------------*/
struct ModelDef  Clients[MMD_MAX_MODEL];

/* ------------------------------------------------------------------------*/
/* ROUTINE DECLARATION                                                     */
/* ------------------------------------------------------------------------*/
void MMDc_P_Init (my_int *ChildId, my_int *model_comm, my_int *inter_comm
		  , my_int *npes);
FCALLSCSUB4(MMDc_P_Init, MMDC_P_INIT, mmdc_p_init, my_PINT, my_PINT, my_PINT
	     , my_PINT);

void MMDc_P_SetInd_and_Mem(my_int *ChildId, my_int *bufsize);
FCALLSCSUB2(MMDc_P_SetInd_and_Mem, MMDC_P_SETIND_AND_MEM, mmdc_p_setind_and_mem
	    , my_PINT, my_PINT);

void MMDc_P_SetInd_and_Mem_2way(my_int *ChildId, my_int *bufsize
				, my_int *parbufsize);
FCALLSCSUB3(MMDc_P_SetInd_and_Mem_2way, MMDC_P_SETIND_AND_MEM_2WAY
	    , mmdc_p_setind_and_mem_2way, my_PINT, my_PINT, my_PINT);

void MMDc_P_GetWaitTime(my_int *ChildId, double *WaitTime);
FCALLSCSUB2(MMDc_P_GetWaitTime, MMDC_P_GETWAITTIME, mmdc_p_getwaittime, my_PINT 
	    ,PDOUBLE);

void MMDc_P_FillBuffer(my_int *ChildId, my_int *PeId, double *array);
FCALLSCSUB3(MMDc_P_FillBuffer, MMDC_P_FILLBUFFER, mmdc_p_fillbuffer, my_PINT 
	    , my_PINT , PDOUBLE);

void MMDc_P_SetBarrier(my_int *ChildId);
FCALLSCSUB1(MMDc_P_SetBarrier, MMDC_P_SETBARRIER, mmdc_p_setbarrier, my_PINT);

void MMDc_P_GetBuffer(my_int *ChildId, my_int *PeId, double *array);
FCALLSCSUB3(MMDc_P_GetBuffer, MMDC_P_GETBUFFER, mmdc_p_getbuffer, my_PINT
	    , my_PINT , PDOUBLE);

void MMDc_P_FreeMem(my_int *ChildId);
FCALLSCSUB1(MMDc_P_FreeMem, MMDC_P_FREEMEM, mmdc_p_freemem, my_PINT);
/* ------------------------------------------------------------------------*/

void MMDc_P_Init (my_int *ChildId, my_int *model_comm, my_int *inter_comm
		  , my_int *npes) {

  struct ModelDef     *ac;     /* actual client */
  struct BufDef       *bPE;

/* convert communicator to C */

  Clients[*ChildId-1].model_comm = MPI_Comm_f2c(*model_comm);
  Clients[*ChildId-1].inter_comm = MPI_Comm_f2c(*inter_comm);

/* get rank and size */
  MPI_Comm_rank (Clients[*ChildId-1].model_comm
		 , &Clients[*ChildId-1].model_rank);
  MPI_Comm_size (Clients[*ChildId-1].model_comm
		 , &Clients[*ChildId-1].model_npes);
  MPI_Comm_remote_size (Clients[*ChildId-1].inter_comm
			, &Clients[*ChildId-1].inter_npes);

/* intra communicater is used for MPI_Get */
  MPI_Intercomm_merge (Clients[*ChildId-1].inter_comm, 0
		       , &Clients[*ChildId-1].intra_comm);

  ac = &Clients[*ChildId-1];

/* allocate structure for client PE information */
 
  *npes = ac->inter_npes;

  bPE = (struct BufDef *) malloc (ac->inter_npes * sizeof(struct BufDef));
  ac->buf     = bPE;

  return;
}

/* ------------------------------------------------------------------------*/

void MMDc_P_SetInd_and_Mem(my_int *ChildId, my_int *bufsize) {
  struct ModelDef     *ac;          /* actual client */
  struct BufDef       *bPE;
  int                 ip;
  int                 tag;

  MPI_Aint            bufsizesum;
  long                index;
  double              *base_ptr;
  MPI_Request         req[1024];
  long                index_buf[1024];
  int                 rCount;

  ac  = &Clients[*ChildId-1];
  bPE = ac->buf;

  bufsizesum = 8;
  index      = 0;
  rCount     = 0;
  tag        = 200;

  /* first stride, compute size and set index */
  for (ip=0;ip<ac->inter_npes; ip++) {    /* loop over all client PEs */
    bPE->BufLen=*bufsize;

    /*  send index to client */
    /*  fprintf(stderr,"SERV BEF SEND  PE %d  INDEX %ip \n", ip, index);*/ 

    index_buf[rCount] = index;    /* preserve value between Isend and Wait */
    MPI_Isend (&index_buf[rCount], 1, MPI_LONG, ip, tag, ac->inter_comm, &req[rCount]);

    /* Isend to avoid hang */
    rCount++;

    bPE->BufIndex = index;

    bufsizesum += bPE->BufLen;
    index      += bPE->BufLen;

    if(rCount == 1024) {        /* maximum of 1024 outstanding requests */
       MPI_Waitall (rCount, req, MPI_STATUSES_IGNORE);
       rCount = 0;
    }

    bPE++;
    bufsize++;
  }

  if(rCount > 0) {              /* synchronize Isend requests           */
    MPI_Waitall (rCount, req, MPI_STATUSES_IGNORE);
  }

/* allocate buffer memory */
  ac->TotalBufferSize = bufsizesum*sizeof(double);

  /*  fprintf (stderr,"MPI_Alloc_mem Parent %d %d\n",ac->model_rank,ac->TotalBufferSize);*/
  MPI_Alloc_mem (ac->TotalBufferSize, MPI_INFO_NULL, &base_ptr);


/* create RMA (one sided communication) window for data buffer */
  
  MPI_Win_create (base_ptr, ac->TotalBufferSize, sizeof(double), MPI_INFO_NULL
		  , ac->intra_comm, &ac->BufWin);
  MPI_Win_fence (0, ac->BufWin);
  MPI_Barrier (ac->intra_comm);

    /* fprintf (stderr,"MPI_Alloc_mem Parent %d %d %x\n",ac->model_rank
       ,ac->TotalBufferSize, base_ptr); */

/* second stride, set buffer pointer */

  bPE = ac->buf;
  for (ip=0;ip<ac->inter_npes; ip++) {     /* loop over all client PEs */

    bPE->SendBuf = base_ptr +(bPE->BufIndex);

    bPE++;
  }

  return;
}

/* ------------------------------------------------------------------------*/

void MMDc_P_SetInd_and_Mem_2way(my_int *ChildId, my_int *bufsize, my_int *parbufsize) {
  struct ModelDef     *ac;          /* actual client */
  struct BufDef       *bPE;
  int                 ip;
  int                 tag;

  MPI_Aint            bufsizesum;
  int                 maxbufsize;
  long                index;
  double              *base_ptr;
  MPI_Request         req[1024];
  long                index_buf[1024];
  int                 rCount;

  ac   = &Clients[*ChildId-1];
  bPE  = ac->buf;

  bufsizesum  = 8;
  index       = 0;
  rCount      = 0;
  tag         = 200;

  /* first stride, compute size and set index */
  for (ip=0;ip<ac->inter_npes; ip++) {    /* loop over all client PEs */
    bPE->BufLen=*bufsize;
    bPE->pBufLen=*parbufsize;

    /*  send index for child coupling to client */
    /*  fprintf(stderr,"SERV BEF SEND  PE %d  INDEX %ip \n", ip, index);*/ 

    index_buf[rCount] = index;    /* preserve value between Isend and Wait */
    MPI_Isend (&index_buf[rCount], 1, MPI_LONG, ip, tag, ac->inter_comm, &req[rCount]);
    /* Isend to avoid hang */
    rCount++;

    bPE->BufIndex = index;
 
    if (*bufsize > *parbufsize){
      maxbufsize= *bufsize;
      /* fprintf (stderr,"maxbuf 1 %d %d BufSize %d %d BufLen %d %d\n",ip,maxbufsize,*bufsize,*parbufsize,bPE->BufLen,bPE->pBufLen );*/
    }else{
      maxbufsize = *parbufsize;
      /* fprintf (stderr,"maxbuf 2 %d %d BufSize %d %d BufLen %d %d\n",ip,maxbufsize,*bufsize,*parbufsize,bPE->BufLen,bPE->pBufLen );*/
      /*  fprintf (stderr,"maxbuf 2 %d %d\n",ip,maxbufsize);*/
    }
    
    bufsizesum += maxbufsize;
    index      += maxbufsize;

    if(rCount == 1024) {        /* maximum of 1024 outstanding requests */
       MPI_Waitall (rCount, req, MPI_STATUSES_IGNORE);
       rCount = 0;
    }
    
    /* fprintf(stderr," SetInd BufLen %d %d, bufInd %d \n", bPE->BufLen
       , bPE->pBufLen, bPE->BufIndex);*/


    bPE++;
    bufsize++;
    parbufsize++;
  }

  if(rCount > 0) {              /* synchronize Isend requests           */
    MPI_Waitall (rCount, req, MPI_STATUSES_IGNORE);
  }

  /* determine buffer size */
    ac->TotalBufferSize = bufsizesum*sizeof(double);

    /*fprintf (stderr,"MPI_Alloc_mem Parent %d %d\n",ac->model_rank,ac->TotalBufferSize);*/
  /* allocate buffer memory */

  MPI_Alloc_mem (ac->TotalBufferSize, MPI_INFO_NULL, &base_ptr);


/* create RMA (one sided communication) window for data buffer */
  
  MPI_Win_create (base_ptr, ac->TotalBufferSize, sizeof(double), MPI_INFO_NULL
		  , ac->intra_comm, &ac->BufWin);
  MPI_Win_fence (0, ac->BufWin);
  MPI_Barrier (ac->intra_comm);

  /*   fprintf (stderr,"MPI_Alloc_mem Parent %d %d %x\n",ac->model_rank
       ,ac->TotalBufferSize, base_ptr); */

/* second stride, set buffer pointer */

  bPE  = ac->buf;
  for (ip=0;ip<ac->inter_npes; ip++) {     /* loop over all client PEs */

    bPE->SendBuf  = base_ptr +(bPE->BufIndex);

    bPE++;
  }

  return;
}

/* ------------------------------------------------------------------------*/

void MMDc_P_GetWaitTime(my_int *ChildId, double *WaitTime){

  struct ModelDef     *ac;                 /* actual client */
  double              t1, t2;

  ac  = &Clients[*ChildId-1];

  *WaitTime = -1.0;

  t1 = MMDc_U_Time();
  /* if(ac->model_rank == 0) fprintf(stderr,"Parent U_TIME_start %d %d %f\n", ac->model_rank, *ChildId, t1); */
  MPI_Barrier (ac->intra_comm);            /* Wait for buffer empty   */
  t2 = MMDc_U_Time();
  /* if(ac->model_rank == 0) fprintf(stderr,"Parent U_TIME end %d %d %f \n", ac->model_rank, *ChildId, t2); */

  *WaitTime = t2-t1;
  return;
}

/* ------------------------------------------------------------------------*/

void MMDc_P_FillBuffer(my_int *ChildId, my_int *PeId, double *array){

  struct ModelDef     *ac;                 /* actual client */
  struct BufDef       *bPE; 
  double              *sntbuf;
  int                 i, ip;

  ac  = &Clients[*ChildId-1];

  ip  = *PeId;
  bPE = &ac->buf[ip];
 
  sntbuf = bPE->SendBuf;
  for (i=0; i<bPE->BufLen; i++){
    *sntbuf++ = array[i];
    /*fprintf(stderr, "%d PARENT FILL ARRAY %d %f %f\n", *PeId, i, array[i], *sntbuf);*/
  }

  bPE++;
  return;
}

/* ------------------------------------------------------------------------*/

void MMDc_P_SetBarrier(my_int *ChildId) {
  
  struct ModelDef     *ac;                     /* actual client */

  ac  = &Clients[*ChildId-1];

  MPI_Barrier (ac->intra_comm);                /* buffer is full */

  return;
}
/* ------------------------------------------------------------------------*/

/* ------------------------------------------------------------------------*/

void MMDc_P_GetBuffer(my_int *ChildId, my_int *PeId, double *array) {
  struct ModelDef     *ac;                 /* actual client */
  struct BufDef       *bPE;
  double              *getbuf;
  int                 nr;
  int                 ip, i;

  ac    = &Clients[*ChildId-1];                 
  ip    = *PeId;
  bPE   = &ac->buf[ip];

  getbuf = bPE->SendBuf;
  nr     = bPE->pBufLen;

  if(nr > 0) {
    for (i=0; i<bPE->pBufLen; i++){
      array[i] = *getbuf;
      /*fprintf(stderr, "%d PARENT ARRAY %d %f %f\n", *PeId+4, i, array[i], *getbuf);*/
      getbuf++;
    }
  }  
  /*fprintf (stderr, " %d after parent C  %d %d %x\n",  *PeId+4, bPE->BufIndex, bPE->pBufLen, bPE->SendBuf);*/
  bPE++;
  return;
}

/* ------------------------------------------------------------------------*/

void MMDc_P_FreeMem(my_int *ChildId){

  struct ModelDef *ac;  /* actual child */
  struct BufDef   *bPE;

  ac = &Clients[*ChildId-1];

  MPI_Win_free(&ac->BufWin);
  MPI_Free_mem(ac->buf->SendBuf);

  free(ac->buf);
 
 return; 
}
/* ------------------------------------------------------------------------*/
