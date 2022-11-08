module structure_generator
    implicit none

    contains
        subroutine generate_sample(seed_value, filename, mean_d, std_d, poro, nx, d_equivalent) 
            ! Generate a 2D fibre structure using random values
            implicit none
            integer :: seed_value
            character :: filename
            real :: mean_d
            real :: std_d
            real :: poro
            integer :: nx
            ! Local definitions
            integer :: i
            integer :: seed_size
            integer, allocatable :: seed(:)
            integer :: nb_fiber
            real :: dist_d
            real :: poro_eff
            real :: d_equivalent

            ! Set the seed
            call random_seed(size=seed_size)
            allocate(seed(seed_size))
            do i=1,seed_size
                seed(i) = seed_value * i
            end do
            call random_seed(put=seed)

            call fibre_distribution(mean_d, std_d, poro, nx, nb_fiber, dist_d, poro_eff, d_equivalent)
        end subroutine generate_sample

        subroutine fibre_distribution(mean_d, std_d, poro, nx, nb_fiber, dist_d, poro_eff, d_equivalent)
            real :: mean_d
            real :: std_d
            real :: poro
            integer :: nx
            integer :: nb_fiber
            real :: dist_d
            real :: poro_eff
            real :: d_equivalent
            integer :: i
            real :: print_value

            do i=1,20
                call norm(mean_d, std_d, print_value)
                print *, print_value
            end do
        end subroutine fibre_distribution

        subroutine norm(mean, std, output)
            ! Generate a sample for normal distribution of mean and std
            real :: mean
            real :: std
            real :: output
            real :: u1, u2
            real, parameter :: pi=3.1415926535

            ! Generate two random numbers
            call random_number(u1)
            call random_number(u2)

            ! Box mueller transform
            output = sqrt(-2*log(u1))*cos(2*pi*u2)
            output = std*output + mean
        end subroutine norm

end module structure_generator

program test
    use structure_generator
    implicit none

    real :: d_eq
    integer :: i
    real :: random_num

    call generate_sample(seed_value=2, filename='file.tiff', mean_d=12.5, std_d=2.85, poro=0.9,  nx=200,  d_equivalent=d_eq)

    do i=1,20
        call random_number(random_num)
        !print *, random_num
    end do

end program test
    
