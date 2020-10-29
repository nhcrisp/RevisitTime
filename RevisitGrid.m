% RevisitGrid is a full-featured example for usage of the RevisitTime
% toolkit that produces a gridded output for different combinations of
% input variables
%
%--- Copyright notice ---%
% Copyright (C) 2017 The University of Manchester
% Written by Nicholas H. Crisp and Sabrina Livadiotti
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
%
%% --- CODE ---%

clear

% Inputs
OneSideView = 0; % Flag for one-side looking [0/1]
DescendAscend = 1; % Flag for access on both descending and ascending pass [0/1]
SSO = 1; % SSO orbit flag

lat = 40; % Target (geocentric) latitude [deg]
dayLimit(1) = 1; % ListPasses Loop Control [days]
dayLimit(2) = 60; % Total Analysis Control [days]

ecc = 0; % Eccentricity
AoP = 0; % Argument of Perigee [deg]

h_min = 700e3; h_max = 800e3; % Min/Max orbit altitude (apogee) [m]
h_delta = 1e3; % Altitude step [m]

if SSO == 0
    inc_min = 60; inc_max = 60; % Min/Max orbit inclination [deg]
    inc_delta = 1; % Inclination step [deg]
end

psi_min = input('Type sensor Min off-nadir FoV (half-cone) angle [deg]:  ');
psi_max = input('Type sensor Max off-nadir FoV (half-cone) angle [deg]:  ');
elv = input('Type elevation angle [deg]:  ');

psi = [psi_min,psi_max];

f_p = ~isempty(psi); 
f_e = ~isempty(elv);   

% Walker constellation configuration
noSats = 6;
noPlanes = 6;
relSpacing = 0;

% Constants
R_E = 6378.137e3; % Earth radius [m]
mu_E = 3.986004418e14; % Earth gravity [m3/s2]
J2 = 1.082629989052e-3; % Earth Zonal Harmonic of Degree 2
omega_E = 7.292115373194000e-05; % Earth Rotation [rad solar second^-1]
Flat = 1/298.257223563; % Earth Flattening Factor
ecc_E = 0.081819221456; % Earth Eccentricity
Tday = 24*3600; % [s]

%% CHECK

if noSats > 1
    if rem(noSats,noPlanes) ~= 0
        warning('Uneven number of satellites per orbital plane')
    end
    
    if noPlanes > noSats
        warning('Number of orbital planes is greater than number of satellites')
    end
end

%% ORBIT PARAMETERS

if f_p==1
    psi_delta = 1; % FoV angle step [deg]
    psi = psi_min:psi_delta:psi_max;
end

% Altitude Vector Determination
h_a = h_min:h_delta:h_max;

% Semimajor axis [m]
a = (R_E + h_a)/(1+ecc); 

% Inclination Vector/Eccentricity/AoP/RAAN Determination for SSO and non-SSO
if SSO == 1 
    dOmegaSS = 360/(365.2421897*Tday)*(pi/180); % SSO Nodal Regression Rate [rad/s]
    ecc = 0; AoP = 0; RAAN = 0;
    inc = nan;
    inc_SS = acosd((-2*a.^(7/2)*dOmegaSS*(1-ecc^2)^2)/(3*R_E^2*J2*sqrt(mu_E))); % Inclination [deg]
else
    inc = inc_min:inc_delta:inc_max;
    RAAN = 0;
    if ecc == 0
        AoP = 0;
    end
    if any(inc < lat)
        cont = input('Inclination may be less than target latitude. Continue? Y/N: ','s');
        if strcmpi(cont, 'N')
            return
        end
    end
end

%% COMPUTATIONS

if f_p == 1
    if SSO == 1 %Sun-Synchronous Orbit
    [X,Y,Z] = ndgrid(psi,a,lat);
    W(:,1) = X(:); W(:,2) = Y(:); W(:,3) = Z(:);
    W(:,4) = repmat(inc_SS',[length(psi),1]);
    elseif SSO == 0 %Non Sun-Synchronous Orbit
    [X,Y,Z,P] = ndgrid(psi,a,lat,inc);
    W(:,1) = X(:); W(:,2) = Y(:); W(:,3) = Z(:); W(:,4) = P(:);
    end
    if f_e == 0
      var = W(:,1);
    elseif f_e == 1
      var = [W(:,1), repmat(elv, [length(W(:,1)),1])];
    end
elseif f_p == 0
    if SSO == 1 %Sun-Synchronous Orbit
    [Y,Z] = ndgrid(a,lat);
    W(:,2) = Y(:); W(:,3) = Z(:); W(:,4) = inc_SS';
    elseif SSO == 0 %Non Sun-Synchronous Orbit
    [Y,Z,P] = ndgrid(a,lat,inc);
    W(:,2) = Y(:); W(:,3) = Z(:); W(:,4) = P(:);
    end
    var = repmat(elv, [size(W,1),1]);
end

% Pre-allocate
maxRevisit = zeros(size(W,1),1);
meanRevisit = zeros(size(W,1),1);
dayLimitOut = zeros(size(W,1),1);

% Loop through inclinations, altitudes, and view angles
wb = waitbar(0,'Please wait...');

for i = 1:size(W,1) 
    coes = [W(i,2), ecc, W(i,4)*(pi/180), RAAN*(pi/180), AoP*(pi/180)];
    [maxRevisit(i),meanRevisit(i), flag] = RevisitCalc(coes,W(i,3),f_p,f_e,noSats,noPlanes,relSpacing,DescendAscend,OneSideView,dayLimit,var(i,:));
    wb = waitbar(i / size(W,1));
end

close(wb)

% Reshaping Max Revisit vector into a length(inc)x length(h_a)x length(psi) matrix 
[MaxRevisit] = ReShape(maxRevisit,flag, SSO, length(inc), length(h_a), length(psi));

maxRevisit_hr = MaxRevisit*24;
maxRevisit_min = MaxRevisit*24*60;
maxRevisit_sec = MaxRevisit*24*60*60;

%% POST - PROCESSING
RevisitPlots
