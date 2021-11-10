module sampler
  implicit none
  integer, parameter :: dp = kind(1.d0)
! This contains all the functions and subroutines called in program
contains

! generates a random r given a
  real(dp) function r_sampler(mode, uniform_rand)
    real(dp) :: mode, uniform_rand, tolerance, r, r_new, diff, func
    integer :: j
    r = 1.d0
    tolerance = 0.d0000000000000001
    do j = 1,100
      diff = 4*(mode**3)*(r**2)*exp(-2*mode*r)
      func = 1 - uniform_rand - exp(-2*mode*r)*(2*(mode**2)*(r**2) + 2*mode*r + 1)
      r_new = r - (func/diff)
      if (abs(r_new - r) < tolerance) then
        r = r_new
        exit
      else
        r = r_new
      end if
    end do
    r_sampler = r
  end function

  real(dp) function phi_sampler(uniform_rand)
    real(dp) :: uniform_rand
    real(dp), parameter :: pi = 4.d0*atan(1.d0)
    phi_sampler = 2.d0 * pi * uniform_rand
  end function

  real(dp) function theta_sampler(uniform_rand)
    real(dp) :: uniform_rand
    theta_sampler = acos(2.d0 * uniform_rand - 1.d0)
  end function

  subroutine histogram(dist,dist_size,resolution,hist)
    integer, intent(in) :: dist_size ! size of input distribution
    real(dp), dimension(dist_size), intent(in) :: dist! input distribution (vector), resolution of histogram
    real(dp), intent(in) :: resolution
    integer, dimension(:), allocatable, intent(out) :: hist ! the histogram itself, a list of integers at each interval
    integer :: j, hist_elem, hist_size
    hist_size = ceiling(maxval(dist)/resolution)
    allocate (hist(hist_size))
    do j = 1,hist_size
      hist(j) = 0
    end do
    do j = 1,dist_size
      hist_elem = ceiling(dist(j)/resolution)
      hist(hist_elem) = hist(hist_elem) + 1
    end do
  end subroutine histogram

end module sampler



program ChargeDensitySampler
  use sampler
  implicit none
  integer :: n_ele, i
  real(dp) :: bohr_rad, atom_num, res
  real(dp), dimension(:,:), allocatable :: ele_dist, ele_dist_rand
  real(dp), dimension(:), allocatable :: r_dist
  integer, dimension(:), allocatable :: r_histogram

! assignments
  bohr_rad = 1.d0
  atom_num = 1.d0
  n_ele = 500000

  allocate (ele_dist(3,n_ele))! arrays use Fortran ordering. (3,n_ele) would be better.
  allocate (ele_dist_rand(3,n_ele))
  call random_number(ele_dist_rand)

! initial electron distribution
  do i = 1,n_ele
    ele_dist(1,i) = r_sampler(atom_num/bohr_rad, ele_dist_rand(1,i))
    ele_dist(2,i) = theta_sampler(ele_dist_rand(2,i))
    ele_dist(3,i) = phi_sampler(ele_dist_rand(3,i))
  end do

! output
  print "(3a20)", "r,", "theta,", "phi "
  do i = 1,n_ele
    print *, i
    print "(3f20.16)", ele_dist(1,i), ele_dist(2,i), ele_dist(3,i)
  end do

  allocate (r_dist(n_ele))
  do i = 1,n_ele
    r_dist(i) = ele_dist(1,i)
  end do
  res = 0.1
  call histogram(r_dist,n_ele,res,r_histogram)

  do i = 1,size(r_histogram)
    print "(i6)", r_histogram(i)
  end do  
  deallocate (ele_dist_rand)
  deallocate (ele_dist)
  deallocate (r_dist)
  deallocate (r_histogram) ! histogram subroutine allocates r_histogram
! make a histogram to compare to expected distribution
! divide r into a grid
end program ChargeDensitySampler
