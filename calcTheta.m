function [ theta, elv_psi, flag, nu ] = calcTheta( lat, coes, DscPass, f_p, f_e, flag, var )
%CALCTHETA Compute theta (half FoV angle) at target latitude
%
% Inputs:
%       lat         Target latitude [deg]
%       coes        Vector of classical orbital elements [sma in m, ecc [-], inc in rad, RAAN in rad, AoP in rad]
%       DscPass     Flag for descending pass viewing
%       f_p         Flag for FoV constraint
%       f_e         Flag for Elv constraint
%       flag        Flag for Elv override
%       var         Variable input vector of psi [deg] and elv [deg]
%
% Outputs:
%       theta       Half FoV angle at target latitude [deg]
%       elv_psi     Elevation of FoV angle
%       flag        Flag for Elv override
%       nu          True anomaly [rad]
%
%--- Copyright notice ---%
% Copyright (C) 2017 The University of Manchester
% Written by Nicholas H. Crisp and Sabrina Livadiotti
%
% This file is part of the RevisitTime toolkit.
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.

%% --- CODE ---%

%To compute Theta and Elv angle in the second case: Psi angle known, Elv angle Unknown

a = coes(1);
ecc = coes(2);
inc = coes(3);
AoP = coes(5);

R_E = 6378.137e3; % Earth radius [m]
Flat = 1/298.257223563; % Earth Flattening Factor
R_Emin = R_E - Flat*R_E;
% Surface radius at target latitude
R_ell = sqrt(((R_E^2*cosd(lat)).^2+(R_Emin^2*sind(lat)).^2)./((R_E*cosd(lat)).^2+(R_Emin*sind(lat)).^2));

% Calculate orbit altitude at target latitude
nu = mod(real(asin(sind(lat)/sin(inc))) - AoP,2*pi); % True anomaly [rad]
u = mod(nu + AoP, 2*pi);
if DscPass == 1
    nu = mod(pi - u - AoP, 2*pi);
end
r_s = (a*(1-ecc^2))./(1 + ecc*cos(nu)); % Orbit radius [m]
h_ell = (R_E - R_ell) + (r_s - R_E); % Orbit altitude [m]

if nargin == 7 && f_p && f_e  % both psi and elv known
    % psi and elv allocation in the respective variables
    psi = var(1,1); 
    elv = var(1,2);
    % Compute theta (half FoV angle) at target latitude, knowing elv
    theta_elv = acosd(cosd(elv)./(1+(h_ell./R_ell)))-elv;
    % Compute theta (half FoV angle) at target latitude, knowing psi
    gam = 180 - asind(((R_ell+h_ell)*sind(psi))./R_ell);
    rho = R_ell*cosd(gam)+(R_ell+h_ell)*cosd(psi); 
    theta_psi = asind((rho*sind(psi))./R_ell); 
    % Compare between theta_elv and theta_psi
    if theta_elv < theta_psi
        theta = theta_elv;
        elv_psi = elv;
        flag = 1; 
    else
        theta = theta_psi;
        elv_psi = psi;
        flag = 0; 
    end
elseif nargin == 7 && f_p && ~f_e %Psi angle is known, Elv angle is unknown
    % psi allocation
    psi = var;
    elv_psi = psi;
    flag = 0;
    % Compute theta (half FoV angle) at target latitude, knowing psi
    gam = 180 - asind(((R_ell+h_ell)*sind(psi))./R_ell);
    rho = R_ell.*cosd(gam)+(R_ell+h_ell).*cosd(psi); %Compute rho (FoV slant range)
    theta = asind((rho*sind(psi))./R_ell); 
elseif nargin == 7 && f_e && ~f_p %Elv angle is known, Psi angle is unknown
    % elv allocation 
    elv = var;
    elv_psi = elv;
    flag = 1;
    % Compute theta (half FoV angle) at target latitude, knowing elv
    theta = acosd(cosd(elv)./(1+(h_ell./R_ell)))-elv;
    %     if deg2rad(lat) > inc
%         theta = 
%     end

elseif nargin == 7 && f_p && f_e 
    % Loop created for theta_range computation in listPasses.m. In case of f_p and f_e = 1, this
    % loop compute theta according to the angle (between psi and elv) that was find out to be more
    % restrictive in the previous loop
    if flag == 1
        elv = var;
        theta = acosd(cosd(elv)./(1+(h_ell./R_ell)))-elv;
        elv_psi = elv;
    elseif flag == 0
        psi = var;
        gam = 180 - asind(((R_ell+h_ell)*sind(psi))./R_ell);
        rho = R_ell*cosd(gam)+(R_ell+h_ell)*cosd(psi); 
        theta = asind((rho*sind(psi))./R_ell);
        elv_psi = psi;
    end
else
    error('Incorrect psi or elv assignment in input arguments')
end
end