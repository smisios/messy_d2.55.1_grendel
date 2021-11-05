MODULE messy_airtraf_tools_ga_parameter

 USE messy_main_constants_mem, ONLY: DP     !HY20140421-3

      IMPLICIT NONE
! ..... General for EA .....
      CHARACTER(LEN=5), PARAMETER      :: version = '1.2.0'
      CHARACTER(LEN=4)                 :: testp
      CHARACTER(LEN=3)                 :: testcode

      INTEGER, PARAMETER               :: imem = 1, istop = 0
      INTEGER, PARAMETER               :: mpop = 200, mdv = 200 ,mobj = 40, mcon = 40, & 
                                          misland = 1, mgend = 200
!HY20140405-2      parameter( iarclimit = mgend + 1 )
      REAL(DP), PARAMETER              :: bv = 1.e-6
      INTEGER, PARAMETER               :: mtpop = mpop*mgend

! ..... Setting
!cnew      parameter( pi = 4.*atan(1.) )
      REAL(DP), PARAMETER              :: pi = 3.14159265358979323846_dp
!      parameter( maxpop = mpop * (ibestn+1) )
      INTEGER, PARAMETER               :: maxpop = mpop*5
      REAL(DP), PARAMETER              :: xminpar = 2.6_dp

!  ...   General for EA   ...
      
      INTEGER                          :: igen, iter, nisland, npop, ndv, nobj, ncon, ndetail, &
                                          imaxmin, igwtobj, nlinear, ngend, ngeocon
      INTEGER                          :: iman, idebug, idman, indveva, indvcnr, ipout
      INTEGER                          :: isetdes, isubar, iadinit, idvsm  !HY20140405-2,idstat 
      REAL(DP)                         :: phub(mdv), phlb(mdv)
      REAL(DP)                         :: geub(misland, mdv), gelb(misland, mdv),&
                                          pave(misland, mdv), stdv(misland, mdv),& 
                                          flat(misland, mdv), pout(misland, mdv) 
      INTEGER                          :: icnrgeo, icnrrnk, icnrvio
      REAL(DP)                 :: gvaltol
      REAL(DP)                 :: objinit, gvalinit
!  ...   Preserve Present Population   ...
      REAL(DP)                 :: xratio
      INTEGER                          :: iarbase
      REAL(DP)                 :: arpopratio, ararcratio
!  ...   Random Number   ...
      INTEGER                          :: idum
      INTEGER                          :: iy
      INTEGER                          :: iv(32)
!  ...   Operator   ...
!      common /share/    sshare, sniche(maxpop)
!HY20140407-1      integer, dimension(:) ::  noffs(misland), nselp(misland),&  
!HY20140407-1                                idvl(misland), isptr(misland)    
      INTEGER                          :: nelite
      INTEGER                          :: nelt(misland) 
      INTEGER                          :: mcros(misland), mmute(misland)
      REAL(DP)                         :: rcros(misland), pcros(misland),&
                                          rmute(misland), pmute(misland),& 
                                          rave(misland), rstd(misland)
      INTEGER                          :: iast(misland), iaint(misland)
      REAL(DP)                         :: rout(misland)
      REAL(DP)                         :: alpdm(mobj,mobj)
      INTEGER                          :: mrpair, loopmax
      REAL(DP)                 :: rnpair
! ..... I/O Data .....
      INTEGER                          :: mtout, mtparam, mtgop, mtgen, mteva, mtobj,&
                                          mtprt, mtall, mtdeb, imtall, irank1, mtsys,&
                                          infopt, imemory
! ..... Parameter .....
      INTEGER                          :: mfiteval, ishare, msdist, iallsh 
      REAL(DP)                 :: shalpha 
      INTEGER                          :: irksp, mprank, ictrlelt, mselection
      REAL(DP)                 :: rkcoef, alpdmpar, arcratio
      INTEGER                          :: ibestn, ialpnorm, iarchiv
      REAL(DP)                 :: arclimit, aratio, clr 
      INTEGER                          :: icrobd, isig2n, ismpop
      INTEGER                          :: mfitsel, ibnarc
! ..... IFLAG .....
      INTEGER                          :: iflag(10)
      INTEGER                          :: istat(misland,20)
      INTEGER                          :: imoga(misland,10)
      INTEGER                          :: nallpop
      REAL(DP)                 :: ptall(mtpop,mdv),gtall(mtpop,mdv),&
                                          fvall(mtpop,mobj),gvall(mtpop,mcon)
      INTEGER                          :: ipinfo(mtpop,4)
      INTEGER                          :: idarc(misland,mtpop),irank(misland,mtpop),&
                                          idparent(misland,maxpop)
      INTEGER                          :: ifitpop(misland,maxpop),&
                                          ifrank(misland,maxpop),nunit(misland,maxpop)
      REAL(DP)                 :: ffit(misland,maxpop),&
                                          sniche(misland,maxpop)
      REAL(DP)                 :: objallmax(misland,mobj),objallmin(misland,mobj),&
                                          objmax(misland,mobj),objmin(misland,mobj),&
                                          objave(misland,mobj),&
                                          gvalallmax(misland,mcon),gvalallmin(misland,mcon),&
                                          gvalmax(misland,mcon),gvalmin(misland,mcon),&
                                          gvalave(misland,mcon)

! ..... OTHERS .....
!      INTEGER :: icheck, ip,ipop, & 
!                 nselect,nfitpop,newpop,numpair,igst,iged,inum,maxl,ioff1,ioff2, &
!                 mate1, mate2,icn,ipst,iped,iofi,iofj,idv,itype,inst,ined, &
!                 nump,nsum1,nsum2,nsum3,iof,idir,new,nmax, &
!                 i1,i2,i3,igs,istobj,iedobj,istcon,iedcon,inew,iarcst,iarced,&
!                 narcin,iarc,irkst,irked,ism,ifp,info,icnr,is,ifinish,nrkp,nrks,&
!                 isub,irk,ntmp,igp,numparent,m1,m2,method,j1,j2,ifix,ioff,ippst,&
!                 ipped,idvst,idved,mt,j,k,ii,iread,id,iall,i,&
!                 iccset,ilnum, ndvt, ncont,nobjt,num,ist,ied,numobj,numcon,numpop 

!      DOUBLE PRECISION :: ptype,gtype,pot,sdv,flo,fup,u,fnor,bd,tmp,f1,g,&
!                          x1,x2,fsum,value,sshare,gs,gj,dini,dijnorm,ftmp,gpop,&
!                          gsub,point,rate,alpha,pmemp1,pmemp2,pm1,pm2,betaq,&
!                          eta,parm,pm,rmu,delta,phlo,phup,g0,s,f2,x,a,gamser,&
!                          gln,gammcf,del,&  
!                          bu,bl,r

END MODULE messy_airtraf_tools_ga_parameter

