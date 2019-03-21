function qd = discord(rho, N_theta, N_phi)
% Requires:     Need TrX.m dm2cm pauliprod2.m
% Author:       Xiao-Ming Lu (luxiaoming@gmail.com)
% Date:         2011/10/22
% License: GPL2
%
% Description: Generate the quantum discord for 4x4 density matrix, 
% the first qubit is the apparatus which is measured. 
%
% Note: all the logarithm functions are with unit log2
%       if you want to find the minimun measurement, use [qd,theta,phi] = discord_a(rho)
%
% Usage: qd = discord(rho <,N_theta <, N_phi >>)
% Parameters:    rho     The density matrix (Require 4by4).
%                N_theta The number of steps for dividing theta.
%                N_phi   The number of steps for dividing phi.
% Return:   qd   the value of quantum discord
%
% Rerefence:  X.-M Lu et. al, Phys. Rev. A 83, 012327 (2011) arXiv:1009.1476  
%  (Calculation of discord is based on the formulas in Appendix B of the paper.)   
%
% measure on the first qubit
% qd = I -cc = VNent(rhoA)+ min_CE - VNent(rho)


%% Copyright (C) Xiao-Ming Lu
%%
%% This program is free software; you can redistribute it and/or
%% modify it under the terms of the GNU General Public License
%% as published by the Free Software Foundation; either version 2
%% of the License, or (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
%% MA 02110-1301, USA.

if nargin < 2
    N_theta = 100;
end
if nargin < 3
    N_phi = N_theta;
end

rhoA = TrX(rho, 2, [2, 2]);
rhoA=rhoA/trace(rhoA);
tau = dm2cm(rho);
a = tau(2:4,1); % 3x1
b = tau(1,2:4)';
R = tau(2:4,2:4);
% the quantity need to minimize is condtion_entropy
conditional_entropy = zeros(N_theta,N_phi); 
    
theta_a = linspace(0,pi,N_theta + 1);
phi_a = linspace(-pi/2,pi/2,N_phi + 1);
theta_a(N_theta + 1) = [];
phi_a(N_phi + 1) = [];

n_array = zeros(3,N_theta,N_phi);
n_array(1,:,:) = sin(theta_a)' * cos(phi_a);
n_array(2,:,:) = sin(theta_a)' * sin(phi_a);
n_array(3,:,:) = cos(theta_a)' * ones(1,N_phi);

f = a'*n_array(:,:);% 1 x N
RTn = R'*n_array(:,:);% 3 x N
gp = sqrt((b(1) + RTn(1,:)) .^2 + (b(2)+RTn(2,:)) .^2 + (b(3) + RTn(3,:)) .^2);
gm = sqrt((b(1) - RTn(1,:)) .^2 + (b(2)-RTn(2,:)) .^2 + (b(3) - RTn(3,:)) .^2);
conditional_entropy(:) =...
            (1+f')/2.*ent_fun(gp'./(1+f')) +...
            (1-f')/2.*ent_fun(gm'./(1-f'));

min_CE = min(min(conditional_entropy));
qd = VNent(rhoA) + min_CE - VNent(rho);

function H = ent_fun(lambda)
	% binary Shannon entropy
	H = -(1 + lambda) ./ 2 .* log((1 + lambda) / 2 + (1 + lambda == 0))...
    	-(1 - lambda) ./ 2 .* log((1 - lambda) / 2 + (1 - lambda == 0));
	H = H / log(2); % chose log(2) as unit

function x = VNent(p)
	% von Neumann 
	e = eig(p);
	x = -e' * log2(e + (e == 0));
