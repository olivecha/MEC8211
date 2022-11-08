module lattice_boltzmann
    implicit none

    contains

        subroutine LBM(poremat, nx, deltaP, dx, d_equivalent, output)
            ! parameter definitions
            integer, dimension(nx, nx) :: poremat
            integer :: nx
            real :: deltaP
            real :: dx
            real :: d_equivalent
            real :: output
            ! Local parameters definitions
            integer :: ny
            real :: omega
            real :: rho0
            real :: mu
            real :: eps
            real :: dt
            integer, dimension(:), allocatable :: solid
            real :: W(9)
            integer :: cx(9)
            integer :: cy(9)
            real, dimension(:, :), allocatable :: N
            real :: FlowRate_old
            real :: FlowRate
            integer :: t_
            integer :: i
            integer :: j

            ! Assign values to parameters
            ny = nx      ! square domain
            omega = 1.0  ! one over relaxation time
            rho0 = 1.0   ! air density
            mu = 1.8e-5  ! air viscosity
            eps = 1e-6   ! convergence criterion

            ! time step
            dt = (1/omega - 0.5) * rho0*dx**2 / 3.0 / mu
            
            ! transforme the fibre mat matrix to 1D vector
            solid = reshape(poremat,[size(poremat)])

            ! weights
            W=(/ 4/9.0,1/9.0,1/36.0,1/9.0,1/36.0,1/9.0,1/36.0,1/9.0,1/36.0 /)
            cx=(/ 0,0,1,1, 1, 0,-1,-1,-1 /)
            cy=(/ 0,1,1,0,-1,-1,-1, 0, 1 /)
            allocate(N(nx*ny, 9))
            do i = 1, nx*ny
                do j = 1, 9
                    N(i,j) = rho0 * W(j)
                end do
            end do

            FlowRate_old = 1.0
            FlowRate = 0.0

            ! temporal loop
            t_ = 1
            ! do while ((abs(FlowRate_old - FlowRate)/FlowRate) >= eps)
                do i=2,9
                    N(:, i) = reshape(cshift(reshape(N(:,i), [nx, ny]), shift=[cx(i), cy(i)], dim=2), [nx*ny])
                end do
                print *, N(300, 1)

        end subroutine LBM


end module

program test
    use structure_generator
    use lattice_boltzmann
    implicit none
    ! Variable declarations
    integer :: seed
    integer :: Nx
    real :: Dx
    real :: deltaP
    real :: mean_D
    real :: std_D
    real :: poro
    real :: d_eq
    character(:), allocatable :: filename
    integer, dimension(:, :), allocatable :: poremat
    real :: output

    ! Variable assignment
    seed     = 2
    Nx       = 200
    Dx       = 1e-6  ! Meters
    deltaP  = 0.1   ! Pa
    mean_D   = 12.5  ! Microns
    std_D    = 2.85  ! Microns
    poro     = 0.9
    filename = 'fiber_mat'
    allocate(poremat(Nx, Nx))

    call generate_sample(seed_value=2, filename=filename, mean_d=mean_D, std_d=std_D, poro=poro,  nx=Nx, &  
                         d_equivalent=d_eq, poremat=poremat)

    call LBM(poremat, Nx, deltaP, Dx, d_eq, output)

end program test




