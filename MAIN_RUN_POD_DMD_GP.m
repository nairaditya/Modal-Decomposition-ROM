clear;clc;close all;
%% POD, DMD and Galerkin Projection model for 2D incompressible flows
% Author: Aditya Nair, University of Nevada Reno
%
% AN provides no guarantees for this code.  Use as-is and for academic
% research use only; no commercial use allowed without permission.  For
% citations, please use the reference below:
%
% Ref 1: A. G. Nair, S. L. Brunton & K. Taira,
% "Networked-oscillator-based modeling and control of unsteady wake flows,"
% Physical Review E, vol 97, 2018
% (https://doi.org/10.1103/PhysRevE.97.063107)
%
% Ref 2: A. G. Nair, K. Taira, B. W. Brunton & S. L. Brunton,
% "Phase-based control of periodic flows,"
% Journal of Fluid Mechanics, vol 927, 2021
% (https://doi.org/10.1017/jfm.2021.735)

% Ref 3: A. G. Nair, S. L. Brunton & K. Taira,
% "Phase-consistent dynamic mode decomposition from multiple overlapping spatial domains,"
% Physical Review Fluids, vol 5, 2020
% (https://doi.org/10.1103/PhysRevFluids.5.074702)
% The code is written for educational clarity and not for speed.

%% Input data
load('Airfoil_aoa9_Re1000_u.mat');
load('Airfoil_aoa9_Re1000_v.mat');
load('Airfoil_aoa9_Re1000_geom.mat');
% 1 fundamental time period of data.
% x_grid - x grid locations
% y_grid - y grid locations
% u - streamwise velocity - stored as 3D tensor (first two columns for
% velocity at grid locations and third for temporal snapshots)
% v - cross stream velocity - stored as 3D tensor (first two columns for
% velocity at grid locations and third for temporal snapshots)
% geom - airfoil geometry
% 1st column - x location of airfoil
% 2nd column - y location of airfoil

% Parameters
Re = 1000;            % Reynolds number
dt = 0.05;            % sampling time between snapshots
m = 360;              % number of grid points in x
n = 360;              % number of grid points in y
Lx = 2;               % length of x-domain
Ly = 2;               % length of y-domain
nt = size(u,3);       % number of snapshots  
t = 0:dt:(nt-1)*0.05; % time vector
nModes = 6;           % number of modes to compute
dx = x_grid(1,2)-x_grid(1,1); % x spacing
dy = y_grid(2,1)-y_grid(1,1); % y spacing

% Compute time-averaged mean
uMean = mean(u,3);vMean = mean(v,3);
u_m = u-uMean;v_m = v-vMean;

% DMD (Dynamic mode decomposition)
UX = reshape(u,[],nt);UY = reshape(v,[],nt);
X = [UX(:,1:nt-1);UY(:,1:nt-1)];
Y = [UX(:,2:nt);UY(:,2:nt)];
[DMD_Modes,Mu,~,b] = DMD(X,Y);
u_DMD = reshape(DMD_Modes(1:(m-1)*(n-1),:),m-1,n-1,[]);
v_DMD = reshape(DMD_Modes((m-1)*(n-1)+1:end,:),m-1,n-1,[]);
% b - dynamic importance
% Mu - eigenvalues of A tilda
% u_DMD - streamwise DMD modes
% v_DMD - crossstream DMD modes

% POD (Proper orthogonal decomposition)
[lambda,a_POD,u_POD,v_POD] = POD(nt,nModes,n-1,m-1,dx,dy,u_m,v_m);
% lamda - energy content
% a_POD - temporal coefficients
% u_POD - streamwise POD modes
% v_POD - crossstream POD modes

% Galerkin projection (compare with a_POD) - for better accuracy, increase
% number of snapshots in one period when performing POD
a_GP = Galerkin_Projection_ROM(x_grid,y_grid,...
    uMean,vMean,u_POD,v_POD,a_POD,n,m,dx,dy,Re,t);

% plot spatial mode
modal_plot(geom, x_grid,y_grid,u_POD(:,:,1));
modal_plot(geom, x_grid,y_grid,u_DMD(:,:,3));

% plot comparison of Galerkin Projection model with POD temporal coeff from
% DNS

figure;
subplot(221);
plot(t,a_POD(:,1),'k-','Linewidth',2);hold on;
plot(t,a_GP(:,1),'r--','Linewidth',2);
xlabel('t');ylabel('a_1');
legend('DNS','GP');
set(gca,'Fontsize',12);
grid on;grid minor;
subplot(223);
plot(t,a_POD(:,3),'k-','Linewidth',2);hold on;
plot(t,a_GP(:,3),'r--','Linewidth',2);
xlabel('t');ylabel('a_3');
legend('DNS','GP');
set(gca,'Fontsize',12);
grid on;grid minor;
