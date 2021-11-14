module sampler
  implicit none
  integer, parameter :: dp = kind(1.d0)

contains

  real(dp) function r_sampler(mode, uniform_rand)
    real(dp) :: mode, uniform_rand, tolerance, r, r_new, diff, func
    integer :: j
    r = 1.d0
    tolerance = 1.d-8
    do j = 1,10000
      diff = 4*(mode**3)*(r**2)*exp(-2*mode*r)
      func = 1.d0 - uniform_rand - exp(-2*mode*r)*(2*(mode**2)*(r**2) + 2*mode*r + 1.d0)
      r_new = r - (func/diff)
      if (abs(r_new - r) < tolerance.and.j>5) then
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

  real(dp) function cos_theta_sampler(uniform_rand) ! Use Cos(theta) for cartesian conversions
    real(dp) :: uniform_rand
    cos_theta_sampler = 2.d0 * uniform_rand - 1.d0
  end function

  subroutine histogram(dist,dist_size,resolution,hist)
    integer, intent(in) :: dist_size ! size of input distribution
    real(dp), dimension(dist_size), intent(in) :: dist! input distribution (vector), resolution of histogram
    real(dp), intent(in) :: resolution
    integer, dimension(:), allocatable :: hist_count ! the histogram itself, a list of integers at each interval
    real(dp), dimension(:), allocatable, intent(out) :: hist
    integer :: j, hist_elem, hist_size
    hist_size = ceiling(maxval(dist)/resolution)
    allocate (hist_count(hist_size))
    hist_count(:) = 0
    do j = 1,dist_size
      hist_elem = ceiling(dist(j)/resolution)
      hist_count(hist_elem) = hist_count(hist_elem) + 1
    end do
    allocate (hist(hist_size))
    do j = 1,hist_size
      hist(j) = dble(hist_count(j)) / dist_size
    end do
    deallocate (hist_count)
  end subroutine histogram

end module sampler


program ChargeDensitySampler
  use sampler
  implicit none
  integer :: n_ele, i
  real(dp) :: bohr_rad, atom_num, res, radius, cos_theta, sin_theta, phi
  real(dp), dimension(:,:), allocatable :: ele_dist, ele_dist_rand, ele_dist_cart
  real(dp), dimension(:), allocatable :: r_dist, r_histogram

  bohr_rad = 1.d0
  atom_num = 1.d0
  n_ele = 500000

  allocate (ele_dist(3,n_ele))! arrays use Fortran ordering. (3,n_ele) would be better.
  allocate (ele_dist_rand(3,n_ele))
  call random_number(ele_dist_rand)

  do i = 1,n_ele
    ele_dist(1,i) = r_sampler(atom_num/bohr_rad, ele_dist_rand(1,i))
    ele_dist(2,i) = acos(cos_theta_sampler(ele_dist_rand(2,i)))
    ele_dist(3,i) = phi_sampler(ele_dist_rand(3,i))
  end do
  deallocate (ele_dist_rand)

! allocate (ele_dist_cart(3,n_ele))
! do i = 1,n_ele
!   radius = r_sampler(atom_num/bohr_rad, ele_dist_rand(1,i))
!   cos_theta = cos_theta_sampler(ele_dist_rand(2,i))
!   sin_theta = sqrt(1 - cos_theta**2)
!   phi = phi_sampler(ele_dist_rand(3,i))
!   ele_dist_cart(1,i) = radius*sin_theta*cos(phi)
!   ele_dist_cart(2,i) = radius*sin_theta*sin(phi)
!   ele_dist_cart(3,i) = radius*cos_theta
! end do

! print "(3a20)", "r,", "theta,", "phi "
! do i = 1,n_ele
!   print *, i
!   print "(3f20.16)", ele_dist_cart(1,i), ele_dist_cart(2,i), ele_dist_cart(3,i)
! end do

  allocate (r_dist(n_ele))
  r_dist(:) = ele_dist(1,:)
  deallocate (ele_dist)
  res = 0.1
  call histogram(r_dist,n_ele,res,r_histogram) ! histogram subroutine allocates r_histogram
  deallocate (r_dist)

  open (file="histogramdata.csv",unit=8,status="replace")
  do i = 1,size(r_histogram)
!   print "(f20.16)", r_histogram(i)
    write (8,"(f20.16,a1,f20.16)") (i-0.5d0)*res,",", r_histogram(i)
  end do  
  close (8)
  deallocate (r_histogram) ! histogram subroutine allocates r_histogram
! deallocate (ele_dist_cart)
end program ChargeDensitySampler
