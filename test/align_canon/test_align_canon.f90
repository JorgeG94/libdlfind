! Standalone regression test for dlf_align_canon_eigvec (src/dlf_coords.f90).
!
! Root cause: dlf_cartesian_align (called from dlf_neb_init to align NEB
! endpoints, src/dlf_neb.f90:278) diagonalises a 3x3 symmetric matrix built
! from the two geometries. For a linear molecule (e.g. the HCN/HNC endpoints
! tested in case 2 below) that matrix is rank 1, so two of its eigenvalues
! are (numerically) degenerate and LAPACK's dsyev is free to return ANY
! orthonormal basis of the corresponding 2-D null space. Which basis comes
! back can change between runs with floating-point noise far below the
! degeneracy -- the motivating report was MKL's dsyev varying under
! concurrent load -- changing the aligned endpoint and hence the whole NEB
! trajectory. dlf_align_canon_eigvec replaces each cluster of near-equal
! eigenvalues with a canonical, deterministic basis built from the
! (numerically stable) projector onto that cluster's subspace, and fixes
! every eigenvector's sign.
!
! This driver exercises that routine directly (case 1: the exact 3x3 matrix
! extracted with gdb from a stalled run; case 2: the same construction
! applied to the actual HCN/HNC endpoint coordinates, converted Angstrom to
! Bohr with CODATA 2018's 0.529177210903, exactly as dlf_cartesian_align
! does internally) under many tiny (1e-15..1e-13) symmetric perturbations,
! plus a well-conditioned bent-water-pair case (case 3) run through the
! full dlf_cartesian_align to confirm the fix leaves a non-degenerate
! alignment unchanged.
!
! On 3c27977 (before the fix) this does not even link: dlf_align_canon_eigvec
! does not exist yet. Comment out CASE1/CASE2's calls into it (replace with
! the raw dlf_matrix_diagonalise output) to see case 1/2 fail outright: the
! "raw" eigenvectors below differ by O(1) between perturbation trials, far
! above the 1e-10 pass threshold used for the canonicalised ones.

program test_align_canon
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob, printl
  implicit none

  integer :: nfail

  nfail = 0

  call case_matrix("case 1: gdb matrix", gdb_matrix(), nfail)
  call case_matrix("case 2: HCN/HNC rsmat", hcn_hnc_rsmat(), nfail)
  call case_bent_water(nfail)

  if (nfail /= 0) then
    write(*,'(A,I0,A)') "FAIL: ", nfail, " assertion(s) failed"
    stop 1
  end if
  write(*,'(A)') "PASS: all assertions ok"

contains

  ! The matrix in the problem report, taken with gdb at the point rsmat is
  ! passed to dlf_matrix_diagonalise inside dlf_cartesian_align for an
  ! HCN/HNC NEB endpoint pair. Eigenvalues 3.18251791, 2.9089e-14,
  ! 2.2109e-16: a near-degenerate pair at the bottom.
  function gdb_matrix() result(m)
    real(rk) :: m(3,3)
    m(1,:) = [ 1.77713427D0, -0.00828267D0, -1.58034390D0]
    m(2,:) = [-0.00828267D0,  0.00003860D0,  0.00736549D0]
    m(3,:) = [-1.58034390D0,  0.00736549D0,  1.40534504D0]
  end function gdb_matrix

  ! Builds the same rsmat = rmat^T * rmat that dlf_cartesian_align builds
  ! internally (src/dlf_coords.f90, dlf_cartesian_align), from the actual
  ! HCN/HNC endpoints under test (test/hcn_minimum.xyz, hnc_minimum.xyz in
  ! the terco repo), converted Angstrom -> Bohr. DL-FIND itself works in
  ! Bohr, so this mirrors what dlf_neb_init actually feeds the aligner.
  function hcn_hnc_rsmat() result(m)
    real(rk) :: m(3,3)
    real(rk), parameter :: ang2bohr = 1.D0/0.529177210903D0
    integer, parameter :: nat = 3
    real(rk) :: coords1(3,nat), coords2(3,nat)
    real(rk) :: center(3), trans(3), rmat(3,3)
    integer :: iat, i, j

    ! HCN, downhill minimum (-distort), test/hcn_minimum.xyz
    coords1(:,1) = [ 1.161096161538D0, -0.001204758176D0, -0.190337130816D0]
    coords1(:,2) = [ 0.347870282512D0, -0.000030399671D0,  0.504910840255D0]
    coords1(:,3) = [-0.528511960915D0,  0.001235157849D0,  1.254153462882D0]
    ! HNC, downhill minimum (+distort), test/hnc_minimum.xyz
    coords2(:,1) = [ 1.131321434060D0,  0.005761457067D0,  1.162840010393D0]
    coords2(:,2) = [-0.559417268123D0, -0.006048726577D0, -0.215517847999D0]
    coords2(:,3) = [ 0.347641351124D0,  0.000287269515D0,  0.523950664661D0]
    coords1 = coords1*ang2bohr
    coords2 = coords2*ang2bohr

    center = 0.D0
    trans  = 0.D0
    do iat = 1, nat
      center = center + coords1(:,iat)
      trans  = trans  + (coords1(:,iat) - coords2(:,iat))
    end do
    center = center/real(nat,rk)
    trans  = trans /real(nat,rk)
    do iat = 1, nat
      coords2(:,iat) = coords2(:,iat) + trans
    end do

    rmat = 0.D0
    do iat = 1, nat
      do i = 1, 3
        do j = 1, 3
          rmat(i,j) = rmat(i,j) + (coords1(i,iat)-center(i))*(coords2(j,iat)-center(j))
        end do
      end do
    end do
    rmat = rmat/real(nat,rk)
    m = matmul(transpose(rmat),rmat)
  end function hcn_hnc_rsmat

  ! Fixed-seed, reproducible "random" symmetric perturbation with entries
  ! of magnitude ~scale (scale itself drawn log-uniformly in [lo,hi] per
  ! trial, as in the task: 1e-15..1e-13). A small self-contained PCG-style
  ! generator is used instead of the intrinsic random_number/random_seed
  ! so the sequence (and hence the test) is identical under gfortran and
  ! nvfortran, whose default seed sizes differ.
  subroutine sym_noise(noise, lo, hi, iseed)
    real(rk), intent(out) :: noise(3,3)
    real(rk), intent(in)  :: lo, hi
    integer,  intent(in)  :: iseed
    integer, parameter :: ik8 = selected_int_kind(18)
    integer(ik8), parameter :: mult = 6364136223846793005_ik8
    integer(ik8), parameter :: incr = 1442695040888963407_ik8
    integer(ik8), parameter :: mask = 9007199254740991_ik8 ! 2**53-1
    integer(ik8) :: state
    real(rk) :: u(7), scale, logscale
    integer :: i, j, n

    state = int(iseed,ik8)*mult + incr
    do n = 1, 7
      state = state*mult + incr
      u(n) = real(iand(ishft(state,-11),mask),rk) / 9007199254740992.D0
    end do

    logscale = log(lo) + u(1)*(log(hi)-log(lo))
    scale = exp(logscale)

    n = 1
    do i = 1, 3
      do j = i, 3
        n = n + 1
        noise(i,j) = (2.D0*u(n) - 1.D0)*scale
        noise(j,i) = noise(i,j)
      end do
    end do
  end subroutine sym_noise

  subroutine case_matrix(name, m0, nfail)
    character(*), intent(in)    :: name
    real(rk),     intent(in)    :: m0(3,3)
    integer,      intent(inout) :: nfail
    integer,  parameter :: ntrials = 300
    real(rk), parameter :: lo = 1.D-15, hi = 1.D-13
    real(rk) :: m(3,3), noise(3,3), eigval(3), raw(3,3), canon(3,3)
    real(rk) :: raw0(3,3), canon0(3,3)
    real(rk) :: diff_raw_max, diff_canon_max
    integer :: t

    m = m0
    call dlf_matrix_diagonalise(3,m,eigval,raw0)
    canon0 = raw0
    call dlf_align_canon_eigvec(3,eigval,canon0)

    diff_raw_max   = 0.D0
    diff_canon_max = 0.D0
    do t = 1, ntrials
      call sym_noise(noise, lo, hi, t)
      m = m0 + noise
      call dlf_matrix_diagonalise(3,m,eigval,raw)
      canon = raw
      call dlf_align_canon_eigvec(3,eigval,canon)
      diff_raw_max   = max(diff_raw_max,   maxval(abs(raw   - raw0)))
      diff_canon_max = max(diff_canon_max, maxval(abs(canon - canon0)))
    end do

    write(*,'(A,A,ES12.4,A,ES12.4)') trim(name), &
      ": max|raw-raw0| = ", diff_raw_max, ", max|canon-canon0| = ", diff_canon_max

    ! Sanity check on the test itself: this case must actually be
    ! ill-conditioned (otherwise the "canon" assertion below would be
    ! trivially satisfied even without the fix, and the test would prove
    ! nothing). If this ever fires, the construction of the matrix case
    ! needs revisiting, not the fix.
    if (diff_raw_max < 1.D-6) then
      write(*,'(A,A)') trim(name), ": test construction error -- " // &
        "raw eigenvectors were not actually unstable; this case does not " // &
        "exercise the bug"
      nfail = nfail + 1
    end if

    ! The actual fix: the canonicalised eigenbasis must be deterministic
    ! and continuous, i.e. stay close to the unperturbed reference for
    ! every perturbation trial.
    if (diff_canon_max > 1.D-10) then
      write(*,'(A,A,ES12.4)') trim(name), &
        ": FAIL -- canonicalised eigenbasis is not stable, max diff = ", diff_canon_max
      nfail = nfail + 1
    end if
  end subroutine case_matrix

  ! Well-conditioned, non-degenerate case: two distinct bent water
  ! geometries (one rotated/translated relative to the other), run through
  ! the full dlf_cartesian_align. The eigenvalues of its internal rsmat are
  ! all well separated here, so dlf_align_canon_eigvec never touches a
  ! cluster of size > 1, and the sign fix alone cannot change the result
  ! (see the "SIGN" note in dlf_align_canon_eigvec's header): the aligned
  ! coordinates must match an independent Kabsch reference to round-off.
  subroutine case_bent_water(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: nat = 3
    real(rk) :: coords1(3,nat), coords2(3,nat), coords2_ref(3,nat)
    real(rk) :: maxdiff
    integer :: iat

    ! Water 1: O at origin-ish, H's forming a bent (104.5 deg) geometry.
    coords1(:,1) = [ 0.000000D0,  0.000000D0,  0.000000D0] ! O
    coords1(:,2) = [ 0.757000D0,  0.586000D0,  0.000000D0] ! H
    coords1(:,3) = [-0.757000D0,  0.586000D0,  0.000000D0] ! H

    ! Water 2: a different bent water, rotated/translated/slightly
    ! distorted relative to water 1, so the two are not already aligned.
    coords2(:,1) = [ 1.200000D0,  2.300000D0,  0.400000D0] ! O
    coords2(:,2) = [ 1.355000D0,  3.034000D0,  0.954000D0] ! H
    coords2(:,3) = [ 0.503000D0,  2.650000D0, -0.312000D0] ! H

    ! Independent SVD-based Kabsch reference for the aligned coordinates,
    ! computed offline with numpy for this exact pair (agrees with a
    ! direct Python transcription of dlf_cartesian_align's own eigenvalue
    ! method to 1e-14, confirming both are solving the same orthogonal
    ! Procrustes problem) and transcribed here.
    coords2_ref(:,1) = [ 0.043382D0, -0.014419D0,  0.000000D0]
    coords2_ref(:,2) = [ 0.764914D0,  0.576418D0,  0.000000D0]
    coords2_ref(:,3) = [-0.808296D0,  0.610000D0,  0.000000D0]

    glob%tatoms = .true.
    if (allocated(glob%spec)) deallocate(glob%spec)
    allocate(glob%spec(nat))
    glob%spec = 1
    printl = 0

    call dlf_cartesian_align(nat, coords1, coords2)

    maxdiff = 0.D0
    do iat = 1, nat
      maxdiff = max(maxdiff, maxval(abs(coords2(:,iat)-coords2_ref(:,iat))))
    end do
    write(*,'(A,ES12.4)') "case 3: bent water pair, max|aligned-ref| = ", maxdiff
    if (maxdiff > 1.D-5) then
      write(*,'(A)') "case 3: FAIL -- aligned coordinates do not match the independent reference"
      nfail = nfail + 1
    end if
  end subroutine case_bent_water

end program test_align_canon
