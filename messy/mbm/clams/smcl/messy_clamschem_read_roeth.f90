Module messy_clamschem_read_roeth

Contains


      subroutine read_roeth(njsr, a, b, c) 
!  This subroutine reads in parameters from Roeth's 3-parameter formula for
!  photolysis rate approximation. 
!
!  last changed 27.04.99, J.U. Grooss

      USE messy_main_constants_mem, ONLY: DP
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: njsr 
      real(kind=DP),dimension(njsr),intent(out)  :: a, b, c
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ios, i 
      character :: dummy*36, line*80
!-----------------------------------------------
!
!     Nullsetzen der Werte
      a=0.
      b=0.
      c=0.

      open(9, file='data/roeth.dat') 
 
!     ueberlesen des Headers (Zeilen mit '*' am Anfang oder nur 0-1 Zeichen)
      read (9, *) line 
      do while(line(1:1) == '*' .or. len_trim(line) <= 1) 
         read (9, '(A)') line
      end do 
      
!     lesen der Werte (erster Wert steht schon auf line)
      
      i = 1 
      read (line, '(A36,3F10.8)',iostat=ios) dummy, a(i), b(i), c(i) 

      do while(ios==0 .and. i<njsr) 
         i = i + 1 
         read (9, '(A36,3F10.8)',iostat=ios) dummy, a(i), b(i), c(i) 
      end do 
 
      close(9) 
      return  
 
      end subroutine read_roeth 

End Module messy_clamschem_read_roeth
