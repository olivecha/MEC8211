%% Sébastien Leclaire (2014) This LBM code was inspired from Iain Haslam (http://exolete.com/lbm/)
%% modified by David Vidal for the purpose of Assignment 3 of the V&V course
%% This code uses the Lattice Boltzmann method (LBM) to compute the flow through the fiber mat.
%% It is to be used as a "black box". Do not try to understand what it does as this is a special implementation of the LBM algorithm
%% and it thus requires an in-deepth knowledge of the LBM to understand it.

function out=LBM(filename,NX,deltaP,dx,d_equivalent)
close all;
NY=NX; % square domain
OMEGA=1.0; % one over relaxation time
rho0=1.0; % density of air
mu=1.8e-5; % viscosity of air

epsilon=1e-6; % convergence criterion for reaching steady state
dt=(1/OMEGA-0.5)*rho0*dx^2/3/mu; % time step

% Reading fiber mat structure in TIFF format
A = imread(filename);
SOLID = reshape(A,1,[]);


W=[4/9,1/9,1/36,1/9,1/36,1/9,1/36,1/9,1/36];
cx=[0,0,1,1, 1, 0,-1,-1,-1];
cy=[0,1,1,0,-1,-1,-1, 0, 1];
N=bsxfun(@times,rho0*ones(NX*NY,9),W);
FlowRate_old=1;
FlowRate=0;

% temporal loop
t_=1;
while (abs(FlowRate_old-FlowRate)/FlowRate)>=epsilon
    for i=2:9
        N(:,i)=reshape(circshift(reshape(N(:,i),NX,NY),[cx(i),cy(i)]),NX*NY,1);
    end
    N_SOLID=N(SOLID,[1 6 7 8 9 2 3 4 5]); % Bounce Back and No Collision
    rho = sum(N,2);
    ux  = sum(bsxfun(@times,N,cx),2)./rho;ux=ux+deltaP/(2*NX*dx*rho0)*dt;
    uy  = sum(bsxfun(@times,N,cy),2)./rho;
    workMatrix=ux*cx+uy*cy;workMatrix=(3+4.5*workMatrix).*workMatrix;
    workMatrix=bsxfun(@minus,workMatrix,1.5*(ux.^2+uy.^2));
    workMatrix=bsxfun(@times,1+workMatrix,W);
    workMatrix=bsxfun(@times,workMatrix,rho);
    N=N+(workMatrix-N)*OMEGA;
    FlowRate_old=FlowRate;
    FlowRate=sum(ux(1:NX))/(NX*dx);
    N(SOLID,:)=N_SOLID;
    t_=t_+1;
end

% Write permeability at the screen, effective porosity and Reynolds number
poro_eff=1-sum(SOLID)/(NX*NY)
Re=rho0*mean(ux(1:NX))*poro_eff*d_equivalent*dx/(mu*(1-poro_eff))
k_in_micron2=mean(ux(1:NX))*mu/deltaP*(NX*dx)/dx^2 %in micron^2

% Vector plot of the flow field 
ux(SOLID)=0; uy(SOLID)=0;ux=reshape(ux,NX,NY)';uy=reshape(uy,NX,NY)';
figure(2);clf;hold on;colormap(gray(2));image(2-reshape(SOLID,NX,NY)');
quiver(1:NX,1:NY,ux,uy,1.5,'b');axis([0.5 NX+0.5 0.5 NY+0.5]);axis image;
title(['Velocity field after ',num2str(t_),' time steps']);
