
#include "mmdc_util.h"

/* ------------------------------------------------------------------------*/
/* VARIABLE DECLARATION                                                    */
/* ------------------------------------------------------------------------*/
struct ModelDef me;

/* ------------------------------------------------------------------------*/
/* ROUTINE DECLARATION                                                     */
/* ------------------------------------------------------------------------*/
void MMDc_C_init (my_int *model_comm, my_int *inter_comm, my_int *npes);
FCALLSCSUB3(MMDc_C_init, MMDC_C_INIT, mmdc_c_init, my_PINT, my_PINT, my_PINT);

void MMDc_C_SetInd_and_Mem(my_int *bufsize);
FCALLSCSUB1(MMDc_C_SetInd_and_Mem, MMDC_C_SETIND_AND_MEM, mmdc_c_setind_and_mem
, my_PINT);

void MMDc_C_SetInd_and_Mem_2way(my_int *bufsize, my_int *parbufsize);
FCALLSCSUB2(MMDc_C_SetInd_and_Mem_2way, MMDC_C_SETIND_AND_MEM_2WAY, mmdc_c_setind_and_mem_2way, my_PINT, my_PINT);

void MMDc_C_GetWaitTime(double *WaitTime);
FCALLSCSUB1(MMDc_C_GetWaitTime, MMDC_C_GETWAITTIME, mmdc_c_getwaittime,PDOUBLE);

void MMDc_C_GetWaitTime_2way(double *WaitTime);
FCALLSCSUB1(MMDc_C_GetWaitTime_2way, MMDC_C_GETWAITTIME_2WAY, mmdc_c_getwaittime_2way,PDOUBLE);

void MMDc_C_GetBuffer(my_int *PeId, double *array);
FCALLSCSUB2(MMDc_C_GetBuffer, MMDC_C_GETBUFFER, mmdc_c_getbuffer 
	    , my_PINT, PDOUBLE);

void MMDc_C_FillBuffer(my_int *PeId, double *array);
FCALLSCSUB2(MMDc_C_FillBuffer, MMDC_C_FILLBUFFER, mmdc_c_fillbuffer 
	    , my_PINT, PDOUBLE);

void MMDc_C_SetBarrier();
FCALLSCSUB0(MMDc_C_SetBarrier, MMDC_C_SETBARRIER, mmdc_c_setbarrier);

void MMDc_C_FreeMem();
FCALLSCSUB0(MMDc_C_FreeMem, MMDC_C_FREEMEM, mmdc_c_freemem);

/* ------------------------------------------------------------------------*/

void MMDc_C_init (my_int *model_comm, my_int *inter_comm, my_int *npes) {

  struct BufDef       *bufPE;

/* convert Communicator to C */

  me.model_comm = MPI_Comm_f2c(*model_comm);
  me.inter_comm = MPI_Comm_f2c(*inter_comm);

/* Get rank and size */
  MPI_Comm_rank (me.model_comm, &me.model_rank);
  MPI_Comm_size (me.model_comm, &me.model_npes);
  MPI_Comm_remote_size (me.inter_comm, &me.inter_npes);

/* Intra communicater is used for MPI_Get */
  MPI_Intercomm_merge (me.inter_comm, 1, &me.intra_comm);

/* Allocate structure for Server PE information */

  bufPE = (struct BufDef *) malloc (me.inter_npes * sizeof(struct BufDef));
  me.buf = bufPE;

  *npes = me.inter_npes; 

  return;
}

/* ------------------------------------------------------------------------*/

void MMDc_C_SetInd_and_Mem(my_int *bufsize) {
  struct ModelDef     *ac;                     /* actual Child */
  struct BufDef       *bPE;
  int                 ip;
   int                 tag;

  MPI_Aint            maxbufsize; 
  long                index;
  double              *base_ptr;

  ac         = &me;
  maxbufsize = 0;
  index      = 0;
    tag      = 200;
/* First stride, Compute size and set index */

  bPE = ac->buf;
  /*      fprintf (stderr,"BEFORE SET MEM IP %d \n", *bufsize); */
  for (ip=0;ip<ac->inter_npes; ip++) {          /* Loop over all server PEs */
    bPE->BufLen=*bufsize;


/*    Receive Index from server */
/*    Preliminary This is only called once, but still may be time-consuming
      on many CPUs with many Arrays */
    /*      tag++;*/
    /*      fprintf (stderr,"BEFORE RECV SET MEM %d %d %d %d\n",bPE->BufLen,ip, tag, *bufsize);*/ 
      MPI_Recv (&index, 1, MPI_LONG, ip, tag, ac->inter_comm
		, MPI_STATUS_IGNORE);

      bPE->BufIndex = index;
      /*      fprintf (stderr,"After  Recv %d %d %d %d \n",ac->model_rank 
	      ,ip,tag,index); */

      if(bPE->BufLen > maxbufsize)  {
         maxbufsize= bPE->BufLen;    /* Determin largest Buffer */
      }

     
    bufsize++;
    bPE++;
  }

  /* Allocate buffer Memory */
  ac->TotalBufferSize = maxbufsize*sizeof(double);       

  /*  fprintf(stderr, "BEFORE ALLOC %d\n",  ac->TotalBufferSize);*/
  /* On Child side, buffer for one pe only */
  MPI_Alloc_mem (ac->TotalBufferSize, MPI_INFO_NULL, &base_ptr);


  /*  fprintf (stderr,"MPI_Alloc_mem Child %d %d %x\n",ac->model_rank
      ,ac->TotalBufferSize, base_ptr); */

/* Create RMA (One Sided Communication) window for data buffer */

  MPI_Win_create (base_ptr, ac->TotalBufferSize, sizeof(double), MPI_INFO_NULL
		  , ac->intra_comm, &ac->BufWin);
  MPI_Win_fence (0, ac->BufWin);
  MPI_Barrier (ac->intra_comm);

  bPE = ac->buf;
  for (ip=0;ip<ac->inter_npes; ip++) {      /* Loop over all child PEs */

      bPE->SendBuf = base_ptr;
      /*fprintf(stderr, "bPE->SendBuf %d\n",bPE->SendBuf);*/

    bPE++;
  }

  return;
}

/* ------------------------------------------------------------------------*/
/* ------------------------------------------------------------------------*/

void MMDc_C_SetInd_and_Mem_2way(my_int *bufsize, my_int *parbufsize) {
  struct ModelDef     *ac;                     /* actual Child */
  struct BufDef       *bPE;
  int                 ip;
  int                 tag;

  MPI_Aint            maxbufsize; 
  MPI_Aint            pmaxbufsize; 
  long                index;
  double              *base_ptr;
  MPI_Request         req[1024];
  long                index_buf[1024];
  int                 rCount;

  ac          = &me;
  maxbufsize  = 0;
  index       = 0;
  tag         = 200;
/* First stride, Compute size and set index */

  bPE  = ac->buf;
  /*bpPE = ac->pbuf;*/
  /*      fprintf (stderr,"BEFORE SET MEM IP %d \n", *bufsize); */
  for (ip=0;ip<ac->inter_npes; ip++) {          /* Loop over all server PEs */
    bPE->BufLen =*bufsize;
    bPE->pBufLen=*parbufsize;


/*    Receive Index from server */
/*    Preliminary This is only called once, but still may be time-consuming
      on many CPUs with many Arrays */
    /*      tag++;*/
    /*      fprintf (stderr,"BEFORE RECV SET MEM %d %d %d %d\n",bPE->BufLen,ip, tag, *bufsize);*/ 
      MPI_Recv (&index, 1, MPI_LONG, ip, tag, ac->inter_comm, MPI_STATUS_IGNORE);

      bPE->BufIndex = index;
      /*      fprintf (stderr,"After  Recv %d %d %d %d \n",ac->model_rank 
	      ,ip,tag,index); */

      if(bPE->BufLen > maxbufsize)  {
         maxbufsize= bPE->BufLen;    /* Determin largest Buffer */
      }
      if(bPE->pBufLen > maxbufsize)  {
         maxbufsize= bPE->pBufLen;    /* Determin largest Buffer */
      }

     if(rCount == 1024) {        /* maximum of 1024 outstanding requests */
       MPI_Waitall (rCount, req, MPI_STATUSES_IGNORE);
       rCount = 0;
    }

    
    bufsize++;
    parbufsize++;
    bPE++;
  }

  /* Allocate buffer Memory */
  ac->TotalBufferSize  = maxbufsize*sizeof(double);       

  /*fprintf(stderr, "BEFORE ALLOC %d %d\n",  ac->model_rank, ac->TotalBufferSize);*/

  /* On Child side, buffer for one pe only */
  MPI_Alloc_mem (ac->TotalBufferSize, MPI_INFO_NULL, &base_ptr);


  /* fprintf (stderr,"MPI_Alloc_mem Child %d %d %x\n",ac->model_rank
     ,ac->TotalBufferSize, base_ptr); */

/* Create RMA (One Sided Communication) window for data buffer */

  MPI_Win_create (base_ptr, ac->TotalBufferSize, sizeof(double), MPI_INFO_NULL
		  , ac->intra_comm, &ac->BufWin);
  MPI_Win_fence (0, ac->BufWin);
  MPI_Barrier (ac->intra_comm);

  bPE = ac->buf;
  /* bpPE = ac->pbuf; */
  for (ip=0;ip<ac->inter_npes; ip++) {      /* Loop over all child PEs */

    bPE->SendBuf  = base_ptr;
    bPE++;
  }

  return;
}

/* ------------------------------------------------------------------------*/

void MMDc_C_GetWaitTime(double *WaitTime) {

  struct ModelDef     *ac;                 /* actual child */
  double              t1,t2;

  ac  = &me;

  t1 = MMDc_U_Time();
  /* if(ac->model_rank == 0) fprintf(stderr,"Child U_TIME_start %d %f \n", ac->model_rank, t1); */

  MPI_Barrier (ac->intra_comm);          /* Wait for server to fill buffer  */
  t2 = MMDc_U_Time();
  /* if(ac->model_rank == 0) fprintf(stderr,"Child U_TIME_end %d %f \n", ac->model_rank, t2); */

  MPI_Barrier (ac->intra_comm);          /* Wait for buffer is filled       */

  *WaitTime = t2-t1;                     /* time waiting for server         */

  return;
}

/* ------------------------------------------------------------------------*/

void MMDc_C_GetWaitTime_2way(double *WaitTime) {

  struct ModelDef     *ac;                 /* actual child */
  double              t1,t2;

  ac  = &me;

  t1 = MMDc_U_Time();
  /* if(ac->model_rank == 0) fprintf(stderr,"Child U_TIME_start %d %f \n", ac->model_rank, t1); */

  MPI_Barrier (ac->intra_comm);          /* Wait for server to fill buffer  */
  t2 = MMDc_U_Time();
  /* if(ac->model_rank == 0) fprintf(stderr,"Child U_TIME_end %d %f \n", ac->model_rank, t2); */

  /*MPI_Barrier (ac->intra_comm);*/          /* Wait for buffer is filled       */

  *WaitTime = t2-t1;                     /* time waiting for server         */

  return;
}

/* ------------------------------------------------------------------------*/

void MMDc_C_FillBuffer(my_int *PeId, double *array){

  struct ModelDef     *ac;                 /* actual child */
  struct BufDef       *bPE; 
  int                 ip, i;
  double              *sntbuf;
  double              *buf;
  int                 nr;

  ac    = &me;                 
  ip    = *PeId;
  bPE  = &ac->buf[ip];

  buf = bPE->SendBuf;
  sntbuf = bPE->SendBuf;
  nr     = bPE->pBufLen;

    if (nr > 0) {

    /*      Get Data from Server */
    for (i=0;i<bPE->pBufLen; i++) {
      /* *sntbuf++ = array[i];*/
      *buf = array[i];
      /*fprintf(stderr, "%d CHILD ARRAY A %d %f %f\n", *PeId, i, array[i], *buf);*/
      buf++;
   }
    
    
    /*  fprintf (stderr,"Before Get %x %d, %d %d %d %d  \n",sntbuf, bPE->SendBuf
	, ac->model_rank, bPE->BufIndex, nr, *PeId);*/

    MPI_Win_lock (MPI_LOCK_SHARED , ip, 0, ac->BufWin);
    MPI_Put (sntbuf, nr, MPI_DOUBLE, ip, bPE->BufIndex, nr, MPI_DOUBLE
	     , ac->BufWin);
    MPI_Win_unlock (ip, ac->BufWin);
    

    }
    /*fprintf (stderr, "%d after schild C %d %d %x \n", *PeId,bPE->BufIndex,bPE->pBufLen,  bPE->SendBuf);*/

  return;
}

/* ------------------------------------------------------------------------*/

void MMDc_C_SetBarrier() {
  
  struct ModelDef     *ac;                     /* actual child */

  ac  = &me;

  MPI_Barrier (ac->intra_comm);                /* buffer is full */

  return;
}

/* ------------------------------------------------------------------------*/
/* ------------------------------------------------------------------------*/

void MMDc_C_GetBuffer(my_int *PeId, double *array) {
  struct ModelDef     *ac;                 /* actual child */
  struct BufDef       *bPE;
  double              *buf;
  int                 nr;
  int                 ip, i;

  ac    = &me;                 
  ip    = *PeId;
  bPE  = &ac->buf[ip];

  buf   = bPE->SendBuf;
  nr    = bPE->BufLen;
  if(nr > 0) {

    /*      Get Data from Server */
    
    /* fprintf (stderr,"Before Get %x %d, %d %d %d %d  \n",buf, bPE->SendBuf
           , ac->model_rank, bPE->BufIndex, nr, *PeId); */

    MPI_Win_lock (MPI_LOCK_SHARED , ip, 0, ac->BufWin);
    MPI_Get (buf, nr, MPI_DOUBLE, ip, bPE->BufIndex, nr, MPI_DOUBLE
	     , ac->BufWin);
    MPI_Win_unlock (ip, ac->BufWin);
    
    /* fprintf (stderr, "after child C %d %d \n", bPE->BufLen, *PeId); */

    for (i=0;i<bPE->BufLen; i++) {
      array[i] = *buf++;
    }
    
  }
  
 return; 
}

/* ------------------------------------------------------------------------*/

void MMDc_C_FreeMem(){

  struct ModelDef *ac;  /* actual child */

  ac    = &me;                 

  MPI_Win_free(&ac->BufWin);
  MPI_Free_mem(ac->buf->SendBuf);

  free(me.buf);
 
 return; 
}
