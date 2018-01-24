      SUBROUTINE foo (a,testmat_in,testmat_in_size,testmat_out,out2)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: a, testmat_in_size
      DOUBLE PRECISION, INTENT(IN) :: testmat_in(testmat_in_size)
      DOUBLE PRECISION, INTENT(OUT) :: testmat_out(testmat_in_size)
      DOUBLE PRECISION, INTENT(OUT) :: out2(5,5)

! Print the input integer
      PRINT*, "Hello from Fortran!"
      PRINT*, "You entered the integer ",a
      PRINT*, ""
      
! Test passing of matrix data
      testmat_out = testmat_in
      testmat_out(1) = 100.
      out2(:,:) = 6.
            
      END SUBROUTINE foo
