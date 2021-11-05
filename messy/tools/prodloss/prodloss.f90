program prodloss

! V. Grewe, DLR-Oberpfaffenhofen, 31 March 2011
! Purpose: Define Identify Prod and Loss-Terms of a family automatically
! 
! Input:  - Definition of Family via namelist
!         - MECCA eqn File: filename via namelist
! Output  - List of equations for production and Loss in a form that it can be used 
!            for diagtrac.tex
! Method: - look at every equation and check whether a member of the family is lost or produced.
!           The net effect is a production or destruction term, respectively.
!           See Crutzen&Schmaizl, 1983
!
! Revision: 26 April 2011: MECCA2 adaptation
implicit none

! parameter definitions
integer, parameter                    ::  nchr=15,nfl=40, nline=250  ! nchr:  Max. Length of character string for chem. species
                                                                     ! nfl:   Max. Length of file names
                                                                     ! nline: Max. length of a line in the eqn-file
integer, parameter                    ::  mfam=30,na=50              ! mfam:  Max. number of family members
                                                                     ! na     Max. number of items per line
integer, parameter                    ::  nerr=100                   ! nerr:  Max. number of erors in eq, file

! Inputs via namelist
character (len=nchr), dimension(mfam) :: fam = ''                    ! Considered family
character (len=nchr)                  :: cfam                        ! Name of family
character (len=nfl)                   :: infile,outfile,infofile     ! Names of input/output files and file for more infos
real, dimension(mfam)                 :: fac                         ! weighting of family members

! Inputs via equation-file
character (len=nline)                 :: lineread                    ! Input line 
character (len=nchr),dimension(na)    :: A                           ! Individual items per line

! Values directly related to input
integer                               :: nfam, nitem,ilast,isel,&    ! nfam: actual number of family members
                                         mecca=2                     ! nitem: actual number of items per line
                                                                     ! ilast: position of last item (=":")  
                                                                     ! isel:  criterium: reation not relevant / prod /loss
                                                                     ! mecca 1 or 2 ?
character (len=nchr)                  :: B                           ! item A(1) adapted to mecca version

! Local
logical                               :: lreact, leof,eol            ! lreact: Equation part in file found
                                                                     ! leof:   End of File
                                                                     ! eol:    End of line
            
integer                               :: i,j                         ! loops
integer                               :: iostat                      ! status of input

! results
real                                  :: sfac,nprod,tfac             ! sfac: +1 for production -1 for loss
                                                                     ! tfac: if actual numbers occur in eqn-file (e.g. 2 NO2)
                                                                     ! nprod: net production
! Input tests
logical, dimension(na)                :: lfound                      ! individual family member occurs in eqn file
logical                               :: allfound                    ! all members occur in eqn file
character (len=nchr), dimension(100)  :: errorineq                   ! error in equ 
integer                               :: ierr                        ! counter for errors
integer                               :: i1,i2 ! op_pj_20130723

namelist /chempl/  cfam,fam,infile,outfile,fac,mecca

! mz_rs_20110412+
! define all values as empty strings:
fam(:) = ""
! mz_rs_20110412-

! read Namelist
READ(*,chempl)

! preset values
lfound(:)=.false.
ierr=0

! Open files
infofile=TRIM(outfile(:nfl-5))//".log"
open(10,file=infile)
! mz_rs_20110412+
! suffix tex added for diagtrac.tex file:
open(20,file=TRIM(outfile)//".tex")
! mz_rs_20110412-open(20,file=outfile)
open(30,file=infofile)

! set values related to input
! Number of Family members
nfam=0
do i=1,mfam
  if (fam(i)>"") nfam=nfam+1
enddo

write(*,*) "          "
write(*,*) "          "
write(*,*) "************************************************"
write(*,*) "* Automatic generation of Prod- and Loss terms *"
write(*,*) "* Input: ",TRIM(infile)
! mz_rs_20110412+
! suffix tex added for diagtrac.tex file:
write(*,*) "* Output: ",TRIM(outfile)//".tex"
! mz_rs_20110412-
write(*,*) "* Diagnostic/Info:",TRIM(infofile)
write(*,*) "*   for MECCA",mecca
write(*,*) "************************************************"

write(*,*) "Family ",TRIM(cfam)," = ",int(fac(1)),"x",TRIM(fam(1)),(" + ",int(fac(i)),"x",TRIM(fam(i)),i=2,nfam)
do j=20,30,10
  write(j,*) "% Automatically generated Prod and loss terms for the family defined below"
  write(j,*) "% Used program: prodloss.f90 (V. Grewe)"
  write(j,*) "% Input: equationfile ",TRIM(infile)
  write(j,*) "% Family: ",TRIM(cfam)," = "
  write(j,*) "%                  ",fac(1)," x ",TRIM(fam(1))
  do i=2,nfam
     write(j,*) "%                + ",fac(i)," x ",TRIM(fam(i))
  enddo
enddo

! Read Eqn-File 
! First step: search for start of equations "#EQUATIONS"
lreact=.false.
do while (.not.lreact) 
 read(10,"(250A)") lineread
 lreact=('#EQUATIONS'==lineread(1:10))
enddo

write(*,*) "... #EQUATION found; start reading equations"

! Second step: read line and disregard, when not reaction
!              delete everything after ":" when line is reaction
leof=.false.
do while (.not.leof)
   A(:)=""
   read(10,FMT="(250A)",IOSTAT=iostat) lineread
   leof=(iostat/=0)
   if (leof) exit
   if (lineread=="") cycle
   if (SCAN(lineread,"=")==0) cycle
   if (lineread(1:1)=="/") cycle
   ilast=SCAN(lineread,":")
   if (ilast==0) cycle
   write(30,*) iostat,lineread
   lineread(ilast+1:)=""
   ! op_pj_20130723+
   ! remove comments like {+ M}; can be more than one!
   DO
      i1 = SCAN(lineread,"{")
      i2 = SCAN(lineread,"}")
      IF ((i1 == 0) .and. (i2 == 0)) exit
      lineread(i1:i2)=""
   END DO
   ! op_pj_20130723-
   i=1
   eol=.false.
! Third Step: Decompose line into items
   do while (.not.eol)
     read(lineread,*,IOSTAT=iostat) A(1:i)
     eol=(iostat/=0)
     i=i+1
     eol=eol.or.(i>na)
   enddo
   nitem=i-2
! Forth Step: Calculate Prod/Loss
   sfac=-1.
   nprod=0.
   tfac=1.
   do i=2,nitem
     select case (A(i)(1:1))
     CASE ("+")
          cycle
     CASE ("=")
          sfac=1.
          cycle
     CASE ("2":"9","0",".","-")
          read(A(i),*) tfac
          cycle
     CASE DEFAULT
       if (SCAN(A(i),"+")/=0.and.(SCAN(A(i),"}")==0)) then
           ierr=ierr+1
           if (ierr>nerr) STOP "Too many errors in equation file"
           errorineq(ierr)=A(1)
       endif
       do j=1,nfam
!       if (fam(j)==A(i)) nprod=nprod+sfac*fac(j)
        if (fam(j)==A(i)) then
          lfound(j)=.true.
          nprod=nprod+sfac*fac(j)*tfac
          ! mz_rs_20110412+
          ! show current contribution, not cumulative value:
          write(30,*) "     -> net-prod= ",sfac*fac(j)*tfac,fam(j),' => total prod = ',nprod
          !write(30,*) "     -> net-prod=",nprod,fam(j)
          ! mz_rs_20110412-write(30,*) "     -> net-prod=",nprod,fam(j)
        endif
       enddo
       tfac=1.
     END SELECT
   enddo  
! Fifth Step: Output
   isel=int(sign(1.,nprod))
   if (nprod==0.) isel=0
   SELECT CASE (isel)
   CASE (0)
      write(30,*) "   -->  not affected"
      write(30,*) "                    "
   CASE (1)
      write(30,*) "   -->  Production term ",nprod
      write(30,*) "                    "
      call adapt_mecca(a(1),mecca,b)
      write(20,*) B," & ", abs(nprod), " Prod"//TRIM(cfam), " & ", (TRIM(A(i))," ",i=2,nitem)
   CASE (-1)
      write(30,*) "   -->  Loss term ",-nprod
      write(30,*) "                    "
      call adapt_mecca(a(1),mecca,b)
      write(20,*) B," & ", abs(nprod), " Loss"//TRIM(cfam)," & ", (TRIM(A(i))," ",i=2,nitem)
   END SELECT
enddo

close(10)
! Sixth Step: Check whether all species are found in equation file to avoid misspellings
allfound=.true. 
do j=1,nfam
   allfound=allfound.and.lfound(j)
   if(.not.lfound(j)) then
       write (*,*) "WARNING: ",fam(j)," NOT FOUND IN EQUATIONS:",TRIM(infile)
       write (30,*) "WARNING: ",fam(j)," NOT FOUND IN EQUATIONS:",TRIM(infile)
   endif
enddo 

if (allfound) then
     write(*,*) "All species found"
     write(30,*) "All species found"
endif
do i=1,ierr
 write(*,*) "Error found in reaction ",errorineq(i)
 write(30,*) "Error found in reaction ",errorineq(i)
enddo

close(20)
write(30,*) "End of File"
close(30)
write(*,*) "End of Program"
write(*,*) "              "


contains
subroutine adapt_mecca(A,m,B)
! Converts the string a into b depending on m
! B=A for mecca1
! replace "{#" and "}" by "<" and ">" for mecca2
! ... and vice versa for MECCA1 ! op_pj_20110505
implicit none
character (len=*), intent(in)     :: A
integer,           intent(in)     :: m
character (len=*), intent(out)    :: B

integer                           :: i=0,j=0

SELECT CASE (m)
CASE (1)
   b=a
! op_pj_20110505+
   i=index(a,'<')
   j=index(a,'>')
   if (i*j.ne.0) b=a(1:i-1)//'{#'//a(i+1:j-1)//'}'//a(j+1:)
!   write(*,*) a,i,j,b
! op_pj_20110505-
   return
CASE(2)
   b=a
   i=index(a,'{')
   j=index(a,'}')
   if (i*j.ne.0) b=a(1:i-1)//'<'//a(i+2:j-1)//'>'//a(j+1:)
!   write(*,*) a,i,j,b ! op_pj_20110505
   return
CASE DEFAULT
 STOP 'Error in Mecca version'
END SELECT

end subroutine adapt_mecca

end program prodloss
