function [ maxRevisit, meanRevisit, flag ] = RevisitCalc( coes, lat, f_p, f_e, noSats, noPlanes, relSpacing, DscPass, OneSideView, dayLimit, var)
%REVISITCALC performs the calculation of revisit time metrics
%
% Inputs:
%       coes            Vector of classical orbital elements [sma in m, ecc [-], inc in rad, RAAN in rad, AoP in rad]        
%       lat             Target latitude [deg]
%       f_p             Flag for FoV constraint
%       f_e             Flag for Elv constraint
%       noSats          Number of satellites
%       noPlanes        Number of equispaced orbital planes
%       relSpacing      Relative spacing of satellites in orbital planes
%       DscPass         Flag for descending pass viewing
%       OneSideView     Flag for one-side viewing (not currently used)
%       dayLimit        Limit for computation (based on memory)
%       var             Variable input vector of FoV [deg] and Elv [deg]
%
% Outputs:
%       maxRevist       Maximum revisit time [days]
%       meanRevisit     Mean revisit time [days]
%       flag            Flag for Elv override
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

flag = [];
reloop = 1;
i = 0;

% Constants
R_E = 6378.137e3; % Earth radius [m]
mu_E = 3.986004418e14; % Earth gravity [m3/s2]
J2 = 1.082629989052e-3; % Earth Zonal Harmonic of Degree 2
omega_E = 7.292115373194000e-05; % Earth Rotation [rad solar second^-1]
Flat = 1/298.257223563; % Earth Flattening Factor
ecc_E = 0.081819221456; % Earth Eccentricity
Tday = 24*3600; % [s]

% Inputs
a = coes(1); % Semi-major Axis [m]
ecc = coes(2); % Eccentricity
inc = coes(3); % Inclination [rad]
RAAN = coes(4); % Right Ascension of Ascending Node [rad]
AoP = coes(5); % Argument of Perigee [rad]

perPlane = noSats/noPlanes;

%% Sensor Properties
lat = atand((1-Flat)^2 * tand(lat)); % Convert geodetic latitude

if lat > inc*(180/pi)
   % warning
end

[theta, elv_psi, flag, nu] = calcTheta(lat, coes, 0, f_p, f_e, flag, var);
% Compute equivalent (full) angular range (longitude) at equator
DlonA = 2*acosd((cosd(theta)-sind(lat)^2)/cosd(lat)^2);

if DscPass == 1
    [theta_2, elv_psi_2, flag_2, nu_2] = calcTheta(lat,coes,1, f_p, f_e, flag, var);
    DlonA_2 = 2*acosd((cosd(theta_2)-sind(lat)^2)/cosd(lat)^2);
end

%% Calculate Orbit
% Calculate orbit period
n = sqrt(mu_E/(a^3)); % Mean Motion [rad s^-1]
p = a*(1-ecc^2); % Semi-parameter [km]
P = 2*pi*sqrt(a^3/mu_E); % Keplerian period [s]
P_n = P * (1/(1 + (3/4)*J2*(R_E/p)^2*(sqrt(1-ecc^2)*(2-3*sin(inc)^2)+(4-5*sin(inc)^2)))); % Nodal period [s]
%P_n = P * (1 - ((3*J2*(4-5*sin(inc)))/(4*(a/R_E)^2*sqrt(1-ecc^2)*(1-ecc*cos(AoP))^2)) - ((3*J2*(1+ecc*cos(AoP))^3)/(2*(a/R_E)^2*(1-ecc^2)^3)));

% Compute drift in longitude crossing at equator due to J2 [rad s^-1]
dRAAN = -(3/2)*n*J2*(R_E/p)^2*cos(inc) + (3/32)*n*J2^2*(R_E^4/p^4)*cos(inc)*(12-4*ecc^2-(80+5*ecc^2)*sin(inc)^2);
Dlon0 = P_n*(-omega_E + dRAAN) * (180/pi); % With J2

%     k_rev = 1:1:((dayLimit*3600*24)/P_n);
%     drev = abs(bsxfun(@minus, abs(k_rev*Dlon0),(2*pi*k_rev)'));
%     [rev_lim,x] = find(drev == min(drev(:)));
%     t_limit = rev_lim * P_n;

% Preallocate array of nodal passes
passNo = 0:1:ceil((dayLimit*3600*24)/P_n);

% Calculate longitude of first pass at ascending node
lon_start = atan2d(cos(AoP+nu)*sin(RAAN)+sin(AoP+nu)*cos(RAAN)*cos(inc), cos(AoP+nu)*cos(RAAN)-sin(AoP+nu)*sin(RAAN)*cos(inc))...
    + (mod(nu+AoP,2*pi)/(2*pi))*P_n*(-omega_E + dRAAN)*(180/pi);

if DscPass == 1
    % Calculate longitude of first pass at descending node
    dLon_2 = atan2d(cos(AoP+nu_2)*sin(RAAN)+sin(AoP+nu_2)*cos(RAAN)*cos(inc), cos(AoP+nu_2)*cos(RAAN)-sin(AoP+nu_2)*sin(RAAN)*cos(inc))...
        + (mod(nu_2+AoP,2*pi)/(2*pi))*P_n*(-omega_E + dRAAN)*(180/pi);
end

while reloop == 1
    i = i+1;
    
    if i == 1
        passNo = 0:1:ceil((dayLimit(1)*3600*24)/P_n);
        passTime_start = P_n*(mod(nu+AoP,360)/(2*pi));
        passTime_0 = passTime_start;
    elseif i > 1
        passTime_start = passTime{i-1}(end)+P_n;
        passNo = passNo(end)+1:1:ceil((dayLimit(1)*i*3600*24)/P_n);
    end
    passTime{i} = passNo*P_n + passTime_0;
    
    lon = mod(lon_start+Dlon0*passNo,360);
    if noPlanes > 1
        Dlon0_mp = [((1:1:noPlanes-1)/noPlanes)'*360; 0];
        dphase = (360*relSpacing)/noSats; % Relative Spacing between planes
        orbfrac_f = ((fliplr(1:1:noPlanes-1).*dphase)/360)';
        Dlon0_f = [orbfrac_f * Dlon0; 0];
        passTimeDiff = [orbfrac_f * P_n; 0];
        lon = mod(repmat(lon, noPlanes, 1) + repmat(Dlon0_mp+Dlon0_f, 1, length(lon)),360);
        passTime{i} = repmat(passTime{i}, noPlanes, 1) + repmat(passTimeDiff, 1, length(passTime{i}));
    end
    if perPlane > 1
        orbfrac_p = ((1:1:perPlane-1)./perPlane)';
        Dlon0_ms = [Dlon0 - orbfrac_p * Dlon0; 0];
        passTimeDiff = [P_n - (orbfrac_p * P_n); 0];
        lon = mod(repmat(lon, perPlane, 1) + repmat(Dlon0_ms, noPlanes, length(lon)),360);
        passTime{i} = repmat(passTime{i}, perPlane, 1) + repmat(passTimeDiff, noPlanes, length(passTime{i}));
    end
    
    if DscPass == 1
        lon_asc = lon;
        lon_2_diff = mod(dLon_2 - lon_start(1,1),360);
        lon_2 = mod(lon + lon_2_diff, 360);
        passTime_2{i} = passTime{i} + P_n*(mod(nu_2-nu,2*pi)/(2*pi));
        lon_all = [lon; lon_2];
    else
        lon_all = lon;
    end
    
    % List contiguous passes by longitude at nodal passing
    all_pass = [0,sort(lon_all(:))',360];
    % Find largest contiguous difference
    maxdiff = max(diff(all_pass));
    
    lon_asc = lon(:)';
    [inview_i{i}, t_view_min_i{i}, t_view_max_i{i}] = listPasses(lat, lon_asc, coes, 0, f_p, f_e, flag, elv_psi, theta, DlonA );
    inview = cell2mat({cat(2, inview_i{:})});
    t_view_min = cell2mat({cat(2, t_view_min_i{:})});
    t_view_max = cell2mat({cat(2, t_view_max_i{:})});
    passTimes = cell2mat({cat(2, passTime{:})});
    
    if DscPass == 1
        lon_dsc = lon_2(:)';
        [inview_i2{i}, t_view_min_i2{i}, t_view_max_i2{i}] = listPasses(lat, lon_dsc, coes, 1, f_p, f_e, flag, elv_psi, theta, DlonA );
        inview_2 = cell2mat({cat(2, inview_i2{:})});
        t_view_min_2 = cell2mat({cat(2, t_view_min_i2{:})});
        t_view_max_2 = cell2mat({cat(2, t_view_max_i2{:})});
        passTimes_2 = cell2mat({cat(2, passTime_2{:})});
        
        inview = [inview, inview_2];
        passTimes = [passTimes, passTimes_2];
        t_view_min = [t_view_min, t_view_min_2];
        t_view_max = [t_view_max, t_view_max_2];
    end
    
    if any(sum(inview,2) < 2)
        reloop = 1;
        if dayLimit(1)*i == dayLimit(2)
            alt = a - R_E; % Altitude [km]
            if flag == 0
                warning(['100% coverage not achieved in time limit (alt = ', num2str(alt/1e3,4), ' km, inc = ', num2str(inc*(180/pi),6), ' deg, FoV = ', num2str(elv_psi,3), ' deg)']);
            elseif flag == 1
                warning(['100% coverage not achieved in time limit (alt = ', num2str(alt/1e3,4), ' km, inc = ', num2str(inc*(180/pi),6), ' deg, ElV = ', num2str(elv_psi,3), ' deg)']);
            end
            maxRevisit = dayLimit(2)+1; meanRevisit= dayLimit(2)+1;
            return
        end
    else
        reloop = 0;
    end
end

% Convert to cell for each row and find non-zero elements
passMatch = cellfun(@find,num2cell(inview,2),'uniformoutput',0);
%
passTimeMax = cellfun(@(x,y) (passTimes(x)+t_view_max(y,x)), passMatch,num2cell([1:1:size(inview,1)]'),'uniformoutput',0);
passTimeMin = cellfun(@(x,y) (passTimes(x)+t_view_min(y,x)), passMatch,num2cell([1:1:size(inview,1)]'),'uniformoutput',0);
pass_max = cellfun(@(x,y) max(diff(sort([x,y]))), passTimeMax, passTimeMin,'uniformoutput',0);
%pass_max = cellfun(@(x) max(diff(sort(passTime(x)))), passMatch,'uniformoutput',0);

pass_avg = cellfun(@(x) mean(diff(sort(passTimes(x)))), passMatch,'uniformoutput',0);
maxRevisit = max(cell2mat(pass_max))/3600/24;
meanRevisit = mean(cell2mat(pass_avg))/3600/24;
end