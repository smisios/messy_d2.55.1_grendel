      SUBROUTINE DGGLSE( M, N, P, A, LDA, B, LDB, C, D, X, WORK, LWORK,
     $                   INFO )
*
*  -- LAPACK driver routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, LWORK, M, N, P
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( * ), D( * ),
     $                   WORK( * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  DGGLSE solves the linear equality-constrained least squares (LSE)
*  problem:
*
*          minimize || c - A*x ||_2   subject to   B*x = d
*
*  where A is an M-by-N matrix, B is a P-by-N matrix, c is a given
*  M-vector, and d is a given P-vector. It is assumed that
*  P <= N <= M+P, and
*
*           rank(B) = P and  rank( ( A ) ) = N.
*                                ( ( B ) )
*
*  These conditions ensure that the LSE problem has a unique solution,
*  which is obtained using a GRQ factorization of the matrices B and A.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrices A and B. N >= 0.
*
*  P       (input) INTEGER
*          The number of rows of the matrix B. 0 <= P <= N <= M+P.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, A is destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)
*          On entry, the P-by-N matrix B.
*          On exit, B is destroyed.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= max(1,P).
*
*  C       (input/output) DOUBLE PRECISION array, dimension (M)
*          On entry, C contains the right hand side vector for the
*          least squares part of the LSE problem.
*          On exit, the residual sum of squares for the solution
*          is given by the sum of squares of elements N-P+1 to M of
*          vector C.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (P)
*          On entry, D contains the right hand side vector for the
*          constrained equation.
*          On exit, D is destroyed.
*
*  X       (output) DOUBLE PRECISION array, dimension (N)
*          On exit, X is the solution of the LSE problem.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= max(1,M+N+P).
*          For optimum performance LWORK >= P+min(M,N)+max(M,N)*NB,
*          where NB is an upper bound for the optimal blocksizes for
*          DGEQRF, SGERQF, DORMQR and SORMRQ.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            LOPT, MN, NR
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY, DGEMV, DGGRQF, DORMQR, DORMRQ,
     $                   DTRMV, DTRSV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      MN = MIN( M, N )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( P.LT.0 .OR. P.GT.N .OR. P.LT.N-M ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, P ) ) THEN
         INFO = -7
      ELSE IF( LWORK.LT.MAX( 1, M+N+P ) ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGGLSE', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Compute the GRQ factorization of matrices B and A:
*
*            B*Q' = (  0  T12 ) P   Z'*A*Q' = ( R11 R12 ) N-P
*                     N-P  P                  (  0  R22 ) M+P-N
*                                               N-P  P
*
*     where T12 and R11 are upper triangular, and Q and Z are
*     orthogonal.
*
      CALL DGGRQF( P, M, N, B, LDB, WORK, A, LDA, WORK( P+1 ),
     $             WORK( P+MN+1 ), LWORK-P-MN, INFO )
      LOPT = WORK( P+MN+1 )
*
*     Update c = Z'*c = ( c1 ) N-P
*                       ( c2 ) M+P-N
*
      CALL DORMQR( 'Left', 'Transpose', M, 1, MN, A, LDA, WORK( P+1 ),
     $             C, MAX( 1, M ), WORK( P+MN+1 ), LWORK-P-MN, INFO )
      LOPT = MAX( LOPT, INT( WORK( P+MN+1 ) ) )
*
*     Solve T12*x2 = d for x2
*
      CALL DTRSV( 'Upper', 'No transpose', 'Non unit', P, B( 1, N-P+1 ),
     $            LDB, D, 1 )
*
*     Update c1
*
      CALL DGEMV( 'No transpose', N-P, P, -ONE, A( 1, N-P+1 ), LDA, D,
     $            1, ONE, C, 1 )
*
*     Sovle R11*x1 = c1 for x1
*
      CALL DTRSV( 'Upper', 'No transpose', 'Non unit', N-P, A, LDA, C,
     $            1 )
*
*     Put the solutions in X
*
      CALL DCOPY( N-P, C, 1, X, 1 )
      CALL DCOPY( P, D, 1, X( N-P+1 ), 1 )
*
*     Compute the residual vector:
*
      IF( M.LT.N ) THEN
         NR = M + P - N
         CALL DGEMV( 'No transpose', NR, N-M, -ONE, A( N-P+1, M+1 ),
     $               LDA, D( NR+1 ), 1, ONE, C( N-P+1 ), 1 )
      ELSE
         NR = P
      END IF
      CALL DTRMV( 'Upper', 'No transpose', 'Non unit', NR,
     $            A( N-P+1, N-P+1 ), LDA, D, 1 )
      CALL DAXPY( NR, -ONE, D, 1, C( N-P+1 ), 1 )
*
*     Backward transformation x = Q'*x
*
      CALL DORMRQ( 'Left', 'Transpose', N, 1, P, B, LDB, WORK( 1 ), X,
     $             N, WORK( P+MN+1 ), LWORK-P-MN, INFO )
      WORK( 1 ) = P + MN + MAX( LOPT, INT( WORK( P+MN+1 ) ) )
*
      RETURN
*
*     End of DGGLSE
*
      END
