module structure_generator
    implicit none

    contains
        subroutine generate_sample(seed_value, filename, mean_d, std_d, poro, nx, d_equivalent, poremat) 
            ! Generate a 2D fibre structure using random values
            implicit none
            integer :: seed_value
            character(len=*) :: filename
            real :: mean_d
            real :: std_d
            real :: poro
            integer :: nx
            ! Local definitions
            integer :: i
            integer :: j
            integer :: k
            integer :: seed_size
            integer, allocatable :: seed(:)
            integer :: nb_fiber
            real, dimension(:), allocatable :: dist_d
            real :: poro_eff
            real :: d_equivalent
            real, dimension(:, :), allocatable :: circle
            integer, dimension(nx, nx) :: poremat
            real, dimension(3) :: c1, c2
            integer :: fiber_count
            integer :: circle_count
            real :: rand_num
            integer :: flag
            real :: di
            real :: xi
            real :: yi
            real :: px
            real :: py

            ! Set the seed
            call random_seed(size=seed_size)
            allocate(seed(seed_size))
            do i=1,seed_size
                seed(i) = seed_value * i
            end do
            call random_seed(put=seed)

            call fibre_distribution(mean_d, std_d, poro, nx, nb_fiber, dist_d, poro_eff, d_equivalent)

            ! redefine the porosity
            poro = poro_eff

            ! define the array storing the circles
            allocate(circle(nb_fiber, 3))
            fiber_count = 1
            circle_count = 1

            ! first circle in the array
            call random_number(rand_num)
            circle(fiber_count, 1) = rand_num * nx
            call random_number(rand_num)
            circle(fiber_count, 2) = rand_num * nx
            circle(fiber_count, 3) = dist_d(fiber_count)

            ! check overlap when creating new circles (fibers)
            do while (fiber_count < nb_fiber)
                
                flag = 0
                di = dist_d(fiber_count+1)
                call random_number(rand_num)
                xi = rand_num * nx
                call random_number(rand_num)
                yi = rand_num * nx
                ! Check overlap 
                do i=1,circle_count
                    c1 = circle(i, 1:3)
                    c2 = (/xi, yi, di/)
                    call check_overlap(c1, c2, nx, flag)
                    if (flag == 1) then
                        exit
                    end if
                end do


                if (flag == 1) then
                    fiber_count = fiber_count + 1
                    continue
                else 
                    fiber_count = fiber_count + 1
                    circle_count = circle_count + 1
                    circle(circle_count, 1) = xi
                    circle(circle_count, 2) = yi
                    circle(circle_count, 3) = di
                end if
            end do 

            ! fill the cells in the porosity grid
            do i = 1, nx
                do j = 1, nx
                    poremat(i, j) = 0
                    px = 0.5 + (i - 1)
                    py = 0.5 + (j - 1)
                    flag = 0
                    do k = 1, circle_count 
                        c1 = circle(k, 1:3)
                        call check_in_circle(px, py, c1, nx, flag)
                        if (flag == 1) then
                            poremat(i, j) = 1
                            exit
                        end if
                    end do
                end do
            end do

           call writeppm1Matrix(poremat, filename)


        end subroutine generate_sample

        subroutine fibre_distribution(mean_d, std_d, poro, nx, nb_fiber, dist_d, poro_eff, d_equivalent)
            implicit none
            real :: pi
            real :: mean_d
            real :: std_d
            real :: poro
            integer :: nx
            integer :: nb_fiber
            real, dimension(:), allocatable :: dist_d
            real :: poro_eff
            real :: poro_eff_old
            real :: d_equivalent
            integer :: i
            real :: dist(10000)
            pi = 3.1415926535

            ! define the fibre normal distribution
            do i=1,10000
                call norm(mean_d, std_d, dist(i))
            end do

            ! define the number of fibres
            nb_fiber = 1

            ! Porosity efficienty
            poro_eff = 1 - sum(dist(1:nb_fiber)**2 /(4*pi))/nx**2

            ! define the number of fibers
            do while (poro_eff >= poro)
                poro_eff_old = poro_eff
                nb_fiber = nb_fiber + 1
                poro_eff = 1 - sum(dist(1:nb_fiber)**2 /(4*pi))/nx**2 
            end do

            if (abs(poro_eff - poro) > abs(poro_eff_old -poro)) then
                nb_fiber = nb_fiber - 1
                poro_eff = poro_eff_old
            end if

            allocate(dist_d(nb_fiber))
            dist_d = dist(1:nb_fiber)
            d_equivalent = sum(dist_d**2) / sum(dist_d)

        end subroutine fibre_distribution

        subroutine norm(mean, std, output)
            ! Generate a sample for normal distribution of mean and std
            implicit none
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

        subroutine check_overlap(c1, c2, nx,  flag)
            implicit none
            ! check if two circles overlap in 9 directions
            integer :: flag
            integer :: nx
            real, dimension(3) :: c1, c2

            if ((c2(1) - c1(1))**2 + (c2(2) - c1(2))**2 < (c2(3) + c1(3))**2) then
                flag = 1
            else if ((c2(1) - c1(1) + nx)**2 + (c2(2) - c1(2))**2 < (c2(3) + c1(3))**2) then
                flag = 1
            else if ((c2(1) - c1(1) - nx)**2 + (c2(2) - c1(2))**2 < (c2(3) + c1(3))**2) then
                flag = 1
            else if ((c2(1) - c1(1))**2 + (c2(2) - c1(1) + nx)**2 < (c2(3) + c1(3))**2) then
                flag = 1
            else if ((c2(1) - c1(1))**2 + (c2(2) - c1(1) - nx)**2 < (c2(3) + c1(3))**2) then
                flag = 1
            else if ((c2(1) - c1(1) + nx)**2 + (c2(2) - c1(2) + nx)**2 < (c2(3) + c1(3))**2) then
                flag = 1
            else if ((c2(1) - c1(1) + nx)**2 + (c2(2) - c1(2) - nx)**2 < (c2(3) + c1(3))**2) then
                flag = 1
            else if ((c2(1) - c1(1) - nx)**2 + (c2(2) - c1(2) + nx)**2 < (c2(3) + c1(3))**2) then
                flag = 1
            else if ((c2(1) - c1(1) - nx)**2 + (c2(2) - c1(2) - nx)**2 < (c2(3) + c1(3))**2) then
                flag = 1
            else
                flag = 0
            end if

        end subroutine check_overlap

        subroutine check_in_circle(px, py, circle, nx,  flag)
            real :: px
            real :: py
            integer :: nx
            real, dimension(3) :: circle
            integer flag
            flag = 0
            if ((px - circle(1))**2 + (py - circle(2))**2 < (circle(3)/2)**2) then
                flag = 1
            else if ((px - (circle(1) + nx))**2 + (py - circle(2))**2 < (circle(3)/2)**2) then
                flag = 1
            else if ((px - (circle(1) - nx))**2 + (py - circle(2))**2 < (circle(3)/2)**2) then
                flag = 1
            else if ((px - circle(1))**2 + (py - (circle(2) - nx))**2  < (circle(3)/2)**2) then
                flag = 1
            else if ((px - circle(1))**2 + (py - (circle(2) + nx))**2  < (circle(3)/2)**2) then
                flag = 1   
            else if ((px - (circle(1) + nx))**2 + ((py - (circle(2) + nx))**2) < (circle(3)/2)**2) then
                flag = 1
            else if ((px - (circle(1) + nx))**2 + ((py - (circle(2) - nx))**2) < (circle(3)/2)**2) then
                flag = 1
            else if ((px - (circle(1) - nx))**2 + ((py - (circle(2) + nx))**2) < (circle(3)/2)**2) then
                flag = 1
            else if ((px - (circle(1) - nx))**2 + ((py - (circle(2) - nx))**2) < (circle(3)/2)**2) then
                flag = 1
            end if
        end subroutine check_in_circle

        subroutine writeppm1Matrix(M,text)
          integer :: M(:,:)
          character(len=*) :: text
          integer :: cols,rows
          integer :: i,j
 
          ! Open File   
          open(unit=100, file=trim(text)//".pbm", status='unknown')
          
          ! Write Header and ppm file type
          write(100,'( A )') "P1"
          write(100,'( A )') "# PPM Type 1 File (generated with fortran)"
 
          ! Write Image Size
          cols = size(M,2)
          rows = size(M,1)
          write(100,'( g0, 1x, g0 )') cols, rows
          
          ! Write Image
          do i=1,rows
            do j=1,cols
              write(100,'( i1 )', advance='no') M(i,j)
            enddo
            write(100,*) ! Endline
          enddo
        end subroutine

end module structure_generator

