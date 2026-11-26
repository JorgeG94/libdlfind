! Standalone regression test for hdlc_make_ut_small (src/dlf_hdlc_hdlclib.f90).
!
! hdlc_make_ut builds the delocalised coordinates' Ut matrix from the
! eigenvectors of G = B Bt, B being the primitive B matrix (nii primitives by
! n Cartesians). When nii > n it now takes them from the n x n matrix Bt B
! instead (hdlc_make_ut_small): u = B v / sqrt(lambda), then re-orthonormalised.
! A cage with highly coordinated centres has tens of thousands of primitives
! for a couple of hundred Cartesians, and the nii x nii route was gigabytes
! and hours there.
!
! The reference here shares no code with either routine: B Bt is formed and
! diagonalised with dsyev in this file, and the projector onto its n6 leading
! eigenvectors is compared with Ut^T Ut. The projector is what has to agree,
! not Ut itself: within a degenerate eigenvalue and in sign the rows are a
! free choice of basis. The rows must also be orthonormal, which is the
! property B v / sqrt(lambda) loses as lambda gets small and the
! Gram-Schmidt passes restore -- case 2 has eigenvalues down to 1e-6, two
! orders above the routine's 1e-8 cutoff, and with the Gram-Schmidt passes
! removed its rows are orthonormal to 8e-12 where the others are to 1e-15,
! which the 1e-12 bound below refuses.
!
! Case 3 calls hdlc_make_ut itself, with more primitives than Cartesians and
! with fewer, to hold the dispatch: both must give the reference's projector.
!
! Build against a libdlfind build tree, e.g.
!   gfortran -I<build>/include test/ut_small/test_ut_small.f90 \
!     -L<build> -ldlfind -llapack -lblas -o test_ut_small

program test_ut_small
  use dlf_parameter_module, only: rk
  use dlf_global, only: printl
  use dlfhdlc_matrixlib, only: matrix, matrix_create, matrix_destroy
  use dlfhdlc_hdlclib, only: hdlc, hdlc_make_ut, hdlc_make_ut_small
  implicit none

  integer :: nfail

  nfail = 0
  printl = 0
  hdlc%internal = .true.
  hdlc%ngroups = 1
  hdlc%contyp = 0

  ! 1: well conditioned, many more primitives than Cartesians
  call one_case('well conditioned, 60 x 24', 60, 24, 1.0_rk, 1.0e-11_rk, .false.)
  ! 2: singular values down to 1e-3, so eigenvalues down to 1e-6
  call one_case('small eigenvalues, 80 x 30', 80, 30, 1.0e-3_rk, 1.0e-7_rk, .false.)
  ! 3: through hdlc_make_ut, on either side of nii = n
  call one_case('dispatch, 500 x 36', 500, 36, 1.0_rk, 1.0e-11_rk, .true.)
  call one_case('dispatch, 20 x 24', 20, 24, 1.0_rk, 1.0e-11_rk, .true.)

  if (nfail == 0) then
    write (*, '(a)') 'test_ut_small: all cases passed'
  else
    write (*, '(a,i0,a)') 'test_ut_small: ', nfail, ' check(s) FAILED'
    stop 1
  end if

contains

  subroutine one_case(label, nii, n, smallest, proj_tol, dispatch)
    ! A B matrix of rank n6 = n - 6 with singular values from 1 down to
    ! `smallest`, through the routine and through the reference
    character(len=*), intent(in) :: label
    integer, intent(in) :: nii, n
    real(rk), intent(in) :: smallest, proj_tol
    logical, intent(in) :: dispatch

    real(rk), parameter :: ORTHO_TOL = 1.0e-12_rk
    type(matrix) :: bprim, ut
    real(rk), allocatable :: b(:, :), left(:, :), right(:, :), g(:, :), w(:), work(:)
    real(rk), allocatable :: pref(:, :), pnew(:, :), gram(:, :)
    real(rk) :: sv, ortho, proj
    integer :: n6, k, i, info, idum

    n6 = n - 6
    if (nii < n6) then
      write (*, '(a)') 'bad case: fewer primitives than degrees of freedom'
      stop 2
    end if

    ! B = L S R^T with orthonormal columns in L (nii x n6) and R (n x n6)
    allocate (left(nii, n6), right(n, n6), b(nii, n))
    call orthonormal_columns(left, 0.37_rk)
    call orthonormal_columns(right, 1.91_rk)
    b = 0.0_rk
    do k = 1, n6
      ! geometric ladder from 1 to `smallest`
      sv = smallest**(real(k - 1, rk)/real(max(n6 - 1, 1), rk))
      do i = 1, n
        b(:, i) = b(:, i) + sv*left(:, k)*right(i, k)
      end do
    end do

    ! reference: the n6 leading eigenvectors of B Bt, by dsyev, here
    allocate (g(nii, nii), w(nii), work(max(1, 34*nii)))
    g = matmul(b, transpose(b))
    call dsyev('V', 'L', nii, g, nii, w, work, size(work), info)
    if (info /= 0) then
      write (*, '(a)') 'dsyev failed in the reference'
      stop 2
    end if
    allocate (pref(nii, nii))
    ! dsyev returns ascending eigenvalues: the leading ones are the last n6
    pref = matmul(g(:, nii - n6 + 1:nii), transpose(g(:, nii - n6 + 1:nii)))

    bprim = matrix_create(nii, n, 'B matrix')
    bprim%data = b
    if (dispatch) then
      call hdlc_make_ut(ut, bprim)
    else
      call hdlc_make_ut_small(ut, bprim, nii, n, n6)
    end if

    allocate (gram(n6, n6), pnew(nii, nii))
    gram = matmul(ut%data, transpose(ut%data))
    do i = 1, n6
      gram(i, i) = gram(i, i) - 1.0_rk
    end do
    ortho = maxval(abs(gram))
    pnew = matmul(transpose(ut%data), ut%data)
    proj = maxval(abs(pnew - pref))

    write (*, '(a,a,2es10.2)') label, ': |Ut Ut^T - 1|, |projector - reference| ', ortho, proj
    if (size(ut%data, 1) /= n6 .or. size(ut%data, 2) /= nii) then
      write (*, '(2x,a)') 'FAILED: Ut has the wrong shape'
      nfail = nfail + 1
    end if
    if (.not. ortho < ORTHO_TOL) then
      write (*, '(2x,a)') 'FAILED: the rows of Ut are not orthonormal'
      nfail = nfail + 1
    end if
    if (.not. proj < proj_tol) then
      write (*, '(2x,a)') 'FAILED: Ut does not span the leading eigenvectors of B Bt'
      nfail = nfail + 1
    end if

    idum = matrix_destroy(bprim)
    idum = matrix_destroy(ut)
  end subroutine one_case

  subroutine orthonormal_columns(a, seed)
    ! Deterministic, well-mixed columns, orthonormalised by Gram-Schmidt twice
    real(rk), intent(out) :: a(:, :)
    real(rk), intent(in) :: seed
    integer :: i, j, k, pass

    do j = 1, size(a, 2)
      do i = 1, size(a, 1)
        a(i, j) = sin(seed + 1.3_rk*real(i, rk)*real(j, rk)) + cos(0.7_rk*real(i, rk) + seed*real(j, rk))
      end do
    end do
    do pass = 1, 2
      do j = 1, size(a, 2)
        do k = 1, j - 1
          a(:, j) = a(:, j) - dot_product(a(:, k), a(:, j))*a(:, k)
        end do
        a(:, j) = a(:, j)/sqrt(dot_product(a(:, j), a(:, j)))
      end do
    end do
  end subroutine orthonormal_columns

end program test_ut_small
