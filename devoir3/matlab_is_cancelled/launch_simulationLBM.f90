! fortran script to launch a fiber structure case
! input variables : 
! 
! seed     : Seed value for the random number generator
! 
! Nx       : Number of grid points on one side of the domain
! 
! Dx       : Grid size
!  
! deltaP   : Pressure drop across porous media
!
! mean_D   : Mean fiber diameter used to generate the structure
!
! std_D    : Standard deviation of the fiber diameters
! 
! poro    : Estimated porosity of the generated fibre structure
! 
! filename : Estimated porosity of the generated fibre structure

program launcher

    implicit none

    ! Variable declarations
    integer :: seed
    integer :: Nx
    real :: Dx
    real :: deltaP
    real :: mean_D
    real :: std_D
    real :: poro
    character(:), allocatable :: filename

    ! Variable assignment
    seed     = 0
    Nx       = 200
    Dx       = 1e-6  ! Meters
    deltaP  = 0.1   ! Pa
    mean_D   = 12.5  ! Microns
    std_D    = 2.85  ! Microns
    poro     = 0.9
    filename = 'fiber_mat.tiff'

    ! Function call
    print *, filename

end program launcher


