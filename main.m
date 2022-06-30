clear

% The construction of the flux functions that we use here is presented in
% detail in A.J. Cerfon and J.P. Freidberg, ``One size fits all" analytic
% solutions to the Grad-Shafranov equation, Physics of Plasmas 17, 032502 (2010)


% Copyright (C) 2010-2012: Antoine Cerfon
% Contact: cerfon@cims.nyu.edu
% 
% This program is free software; you can redistribute it and/or modify 
% it under the terms of the GNU General Public License as published by 
% the Free Software Foundation; either version 2 of the License, or 
% (at your option) any later version.  This program is distributed in 
% the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
% PARTICULAR PURPOSE.  See the GNU General Public License for more 
% details. You should have received a copy of the GNU General Public 
% License along with this program; 
% if not, see <http://www.gnu.org/licenses/>.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Specify equilibrium parameters of interest (in this example, we took
%   an FRC equilibrium as an example)
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

method = 'frc';%Choose experimental device

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = 0;%A as defined in the article, associated with beta regime of interest

switch lower(method)

    case 'tftr'
        
%%%%%%%% TFTR PARAMETERS %%%%%%%%%%%%%%
R0 = 2.5;
epsilon = 0.87/R0;
kappa = 1;
delta  = 0;
D = [0,-0.005,-0.01,-0.015,-0.02,-0.024,-0.027,-0.029,-0.0305];% Contour values for good rendering for A=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'jet'
%%%%%%%% JET PARAMETERS %%%%%%%%%%%%%%
epsilon = 1/3;
kappa = 1.7;
delta  = 0.25;
D = [0,-0.008,-0.016,-0.024,-0.031,-0.036,-0.039,-0.04025];% Contour values for good rendering for A=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'iter'
%%%%%%%% ITER PARAMETERS %%%%%%%%%%%%%%
epsilon = 2/6.2;
kappa = 1.7;
delta  = 0.33;
D = [0,-0.008,-0.016,-0.024,-0.031,-0.035,-0.037];% Contour values for good rendering for A=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'mast'
%%%%%%%% MAST PARAMETERS %%%%%%%%%%%%%%
R0 = 0.85;
epsilon = 0.65/R0;
kappa = 2.45;
delta  = 0.5;
D = [0,-0.025,-0.06,-0.1,-0.145,-0.185,-0.22,-0.238];% Contour values for good rendering for A=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'nstx'
%%%%%%%% NSTX PARAMETERS %%%%%%%%%%%%%%
epsilon = 0.67/0.86;
kappa = 2;
delta  = 0.35;
D = [0,-0.025,-0.06,-0.1,-0.145,-0.185,-0.22,-0.238];% Contour values for good rendering for A=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'frc'
%%%%%%%% FRC PARAMETERS %%%%%%%%%%%%%%
epsilon = 0.99;
kappa = 10;
delta  = 0.7;
D = [0,-0.015,-0.06,-0.14,-0.25,-0.35,-0.43];% Contour values for good rendering for A=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'spheromak';
%%%%%%%% SPHEROMAK PARAMETERS %%%%%%%%%%%%%%
epsilon = 0.95;
kappa = 1;
delta  = 0.2;
D = [0,-0.01,-0.03,-0.06,-0.1,-0.145,-0.185,-0.22,-0.24,-0.252];% Contour values for good rendering for A=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

alpha = asin(delta);%Define alpha as defined in paper in terms of the triangularity delta
curv1 = -(1+alpha)^2/(epsilon*kappa^2);%Curvature at the outer equatorial point
curv2 = -kappa/(epsilon*(cos(alpha))^2);%Curvature at the top
curv3 = (1-alpha)^2/(epsilon*kappa^2);%Curvature at the inner equatorial point

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   Construct linear system associated with boundary conditions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%% Homogeneous solutions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U = [psi1(1+epsilon,0) psi2(1+epsilon,0) psi3(1+epsilon,0) psi4(1+epsilon,0) psi5(1+epsilon,0) psi6(1+epsilon,0) psi7(1+epsilon,0)
    psi1(1-epsilon,0) psi2(1-epsilon,0) psi3(1-epsilon,0) psi4(1-epsilon,0) psi5(1-epsilon,0) psi6(1-epsilon,0) psi7(1-epsilon,0)
    psi1(1-epsilon*delta,kappa*epsilon) psi2(1-epsilon*delta,kappa*epsilon) psi3(1-epsilon*delta,kappa*epsilon) psi4(1-epsilon*delta,kappa*epsilon) psi5(1-epsilon*delta,kappa*epsilon) psi6(1-epsilon*delta,kappa*epsilon) psi7(1-epsilon*delta,kappa*epsilon)
    psi1x(1-epsilon*delta,kappa*epsilon) psi2x(1-epsilon*delta,kappa*epsilon) psi3x(1-epsilon*delta,kappa*epsilon) psi4x(1-epsilon*delta,kappa*epsilon) psi5x(1-epsilon*delta,kappa*epsilon) psi6x(1-epsilon*delta,kappa*epsilon) psi7x(1-epsilon*delta,kappa*epsilon)
    curv1*psi1x(1+epsilon,0)+psi1yy(1+epsilon,0) curv1*psi2x(1+epsilon,0)+psi2yy(1+epsilon,0) curv1*psi3x(1+epsilon,0)+psi3yy(1+epsilon,0) curv1*psi4x(1+epsilon,0)+psi4yy(1+epsilon,0) curv1*psi5x(1+epsilon,0)+psi5yy(1+epsilon,0) curv1*psi6x(1+epsilon,0)+psi6yy(1+epsilon,0) curv1*psi7x(1+epsilon,0)+psi7yy(1+epsilon,0)
    curv3*psi1x(1-epsilon,0)+psi1yy(1-epsilon,0) curv3*psi2x(1-epsilon,0)+psi2yy(1-epsilon,0) curv3*psi3x(1-epsilon,0)+psi3yy(1-epsilon,0) curv3*psi4x(1-epsilon,0)+psi4yy(1-epsilon,0) curv3*psi5x(1-epsilon,0)+psi5yy(1-epsilon,0) curv3*psi6x(1-epsilon,0)+psi6yy(1-epsilon,0) curv3*psi7x(1-epsilon,0)+psi7yy(1-epsilon,0)
    curv2*psi1y(1-epsilon*delta,kappa*epsilon)+psi1xx(1-epsilon*delta,kappa*epsilon) curv2*psi2y(1-epsilon*delta,kappa*epsilon)+psi2xx(1-epsilon*delta,kappa*epsilon) curv2*psi3y(1-epsilon*delta,kappa*epsilon)+psi3xx(1-epsilon*delta,kappa*epsilon) curv2*psi4y(1-epsilon*delta,kappa*epsilon)+psi4xx(1-epsilon*delta,kappa*epsilon) curv2*psi5y(1-epsilon*delta,kappa*epsilon)+psi5xx(1-epsilon*delta,kappa*epsilon) curv2*psi6y(1-epsilon*delta,kappa*epsilon)+psi6xx(1-epsilon*delta,kappa*epsilon) curv2*psi7y(1-epsilon*delta,kappa*epsilon)+psi7xx(1-epsilon*delta,kappa*epsilon)];


%%%%%%%%%%%%% Particular solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B = -[A*psipart1(1+epsilon,0)+(1-A)*psipart2(1+epsilon,0)
    A*psipart1(1-epsilon,0)+(1-A)*psipart2(1-epsilon,0)
    A*psipart1(1-epsilon*delta,kappa*epsilon)+(1-A)*psipart2(1-epsilon*delta,kappa*epsilon)
    A*psipart1x(1-epsilon*delta,kappa*epsilon)+(1-A)*psipart2x(1-epsilon*delta,kappa*epsilon)
    A*(curv1*psipart1x(1+epsilon,0)+psipart1yy(1+epsilon,0))+(1-A)*(curv1*psipart2x(1+epsilon,0)+psipart2yy(1+epsilon,0))
    A*(curv3*psipart1x(1-epsilon,0)+psipart1yy(1-epsilon,0))+(1-A)*(curv3*psipart2x(1-epsilon,0)+psipart2yy(1-epsilon,0))
    A*(curv2*psipart1y(1-epsilon*delta,kappa*epsilon)+psipart1xx(1-epsilon*delta,kappa*epsilon))+(1-A)*(curv2*psipart2y(1-epsilon*delta,kappa*epsilon)+psipart2xx(1-epsilon*delta,kappa*epsilon))];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           Solve linear system for free coefficients c_i
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C = U\B;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           Plot the solution: the flux contours are given by Z
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y] = meshgrid(0:.01:1+epsilon+0.5,-kappa*epsilon-0.1:.01:kappa*epsilon+0.1);
Z = psitot(X,Y,C(1),C(2),C(3),C(4),C(5),C(6),C(7),A);%

figure(1)
contour(X, Y, Z, D,'LineWidth',1.5);
set(gca,'FontName','Times','FontSize',16)
axis equal
xlim([0 1+epsilon+0.5])
ylim([-kappa*epsilon-0.2 kappa*epsilon+0.2])
title(['\epsilon =',num2str(epsilon),', ','\kappa =',num2str(kappa),', ','\delta =',num2str(delta)],'FontName','Times','FontSize',16)
colormap jet