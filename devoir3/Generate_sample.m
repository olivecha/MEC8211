
%% Create a 2D fiber structure and export it in tiff format to be used in the LBM code
function [d_equivalent]= Generate_sample(seed,filename,mean_d,std_d,poro,nx)
%
%INPUT VARIABLES:
%
% SEED: integer representing the seed for initializing the random
% generator. If seed=0, automatic seed generation. If you want to reproduce
% the same fiber structure, use the same seed (fibers will be located at the same place).
%
% MEAN_D: contains the mean fiber to be used
%
% STD_D: contains the standard deviation of the fiber diameters
%
% PORO: estimated porosity of the fiber structure to be generated
%
% NX: domain lateral size in grid cell
%
% OUTPUT VARIABLE:
%
% D_EQUIVALENT: equivalent diameter to be used to represent the fiber size distribution

% initialize seed for random generator
if (seed==0)
    rng('shuffle');  % random seed
else
    rng(seed);  % for reproducibility
end

% Determine distribution of fibers
[nb_fiber,dist_d,poro_eff,d_equivalent]=distribution_of_fiber(mean_d,std_d,poro,nx);

poro=poro_eff;


circle = ones(nb_fiber,3);   %declaring an array to store the shape data (x,y,r)

poremat = zeros(nx,nx);

% positioning fibers
fiber_count = 1; %counter
circle(fiber_count,1) = rand()*nx;   %x-coordinate of centre
circle(fiber_count,2) = rand()*nx;   %y-coordinate of centre
circle(fiber_count,3) = dist_d(fiber_count);  %fiber diameter

%Checking for overlapping with previous fiber (circle)
while (fiber_count < nb_fiber)

    flag = 0;
    di = dist_d(fiber_count);
    xi = rand()*nx;
    yi = rand()*nx;
    %Overlap check in all possible 9 directions for periodicity reason
    for i = 1:fiber_count
        if  (xi - circle(i,1))^2 + (yi - circle(i,2))^2 < (di + circle(i,3))^2 ||...
                (xi - circle(i,1) + nx)^2 + (yi - circle(i,2))^2 < (di + circle(i,3))^2 ||...
                (xi - circle(i,1) - nx)^2 + (yi - circle(i,2))^2 < (di + circle(i,3))^2 ||...
                (xi - circle(i,1))^2 + (yi - circle(i,2) + nx)^2 < (di + circle(i,3))^2 ||...
                (xi - circle(i,1))^2 + (yi - circle(i,2) - nx)^2 < (di + circle(i,3))^2 ||...
                (xi - circle(i,1) + nx)^2 + (yi - circle(i,2) + nx)^2 < (di + circle(i,3))^2 ||...
                (xi - circle(i,1) + nx)^2 + (yi - circle(i,2) - nx)^2 < (di + circle(i,3))^2 ||...
                (xi - circle(i,1) - nx)^2 + (yi - circle(i,2) + nx)^2 < (di + circle(i,3))^2 ||...
                (xi - circle(i,1) - nx)^2 + (yi - circle(i,2) - nx)^2 < (di + circle(i,3))^2

            flag = 1;
            break;
        end
    end

    if(flag == 1)
        continue
    end

    fiber_count = fiber_count + 1;
    circle(fiber_count,1) = xi;
    circle(fiber_count,2) = yi;
    circle(fiber_count,3) = di;
end

%filling in the cell grids
for i = 1:nx
    for j = 1:nx
        px = 0.5 + (i - 1);
        py = 0.5 + (j - 1);
        for k = 1:nb_fiber
            %checking if a grid cell belongs to a circle or its periodic image
            if (px - circle(k,1))^2 + (py - circle(k,2))^2 < (circle(k,3)/2)^2 || ...
                    (px - (circle(k,1) + nx))^2 + (py - circle(k,2))^2 < (circle(k,3)/2)^2 ||...
                    (px - (circle(k,1) - nx))^2 + (py - circle(k,2))^2 < (circle(k,3)/2)^2 ||...
                    (px - circle(k,1))^2 + (py - (circle(k,2) - nx))^2  < (circle(k,3)/2)^2 ||...
                    (px - circle(k,1))^2 + (py - (circle(k,2) + nx))^2  < (circle(k,3)/2)^2 ||...
                    (px - (circle(k,1) + nx))^2 + ((py - (circle(k,2) + nx))^2) < (circle(k,3)/2)^2 ||...
                    (px - (circle(k,1) + nx))^2 + ((py - (circle(k,2) - nx))^2) < (circle(k,3)/2)^2 ||...
                    (px - (circle(k,1) - nx))^2 + ((py - (circle(k,2) + nx))^2) < (circle(k,3)/2)^2 ||...
                    (px - (circle(k,1) - nx))^2 + ((py - (circle(k,2) - nx))^2) < (circle(k,3)/2)^2
                %poremat(i,j) = 0;
                poremat(i,j) = 1;

                break
            end
        end
    end
end

imwrite(logical(poremat),filename,'tiff');

data=imread(filename);
imshow(data);

end

function [nb_fiber,dist_d,poro_eff,d_equivalent]=distribution_of_fiber(mean_d,std_d,poro,nx)

dist=normrnd(mean_d,std_d,[1,10000]);
nb_fiber=1;
poro_eff=1-sum(dist(1:nb_fiber).^2/4*pi)/nx^2;

while poro_eff >= poro
    poro_eff_old=poro_eff;
    nb_fiber=nb_fiber+1;
    poro_eff=1-sum(dist(1:nb_fiber).^2/4*pi)/nx^2;
end
if (abs(poro_eff-poro)>abs(poro_eff_old-poro))
    nb_fiber=nb_fiber-1;
    poro_eff=poro_eff_old;
end

dist_d=dist(1:nb_fiber);

d_equivalent=(sum(dist_d.^2)/sum(dist_d))

end
