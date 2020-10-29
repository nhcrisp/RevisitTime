function  [ inView, t_view_min, t_view_max ] = listPasses( lat, passLon, coes, DscPass, f_p, f_e, flag, elv_psi, theta, DlonA )
%LISTPASSES 
%
% Inputs:
%       lat         Target latitude [deg]
%       pass_lon    Vector of longitudes [deg]
%       coes        Vector of classical orbital elements [sma in m, ecc [-], inc in rad, RAAN in rad, AoP in rad]
%       DscPass     Flag for descending pass viewing
%       f_p         Flag for FoV constraint
%       f_e         Flag for Elv constraint
%       flag        Flag for Elv override   
%       elv_psi     Elevation of FoV angle
%       theta       Half FoV angle at target latitude [deg]
%       DlonA       Equivalent angular range at equator [deg]
%
% Outputs:
%       inView      Index matrix of views
%       t_view_min  Minimum time of view [s]
%       t_view_max  Maximum time of view [s]
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
%
%% --- CODE ---%
R_E = 6378.137e3; % Earth radius [m]
J2 = 1.082629989052e-3; % Earth Zonal Harmonic of Degree 2
mu_E = 3.986004418e14; % Earth gravity [m3/s2]
omega_E = 7.292115373194000e-05; % Earth Rotation [rad solar second^-1]
Flat = 1/298.257223563; % Earth Flattening Factor

a = coes(1); % Semi-major Axis [m]
ecc = coes(2); % Eccentricity
inc_r = coes(3); % Inclination [rad]
RAAN = coes(4); % Right Ascension of Ascending Node
AoP = coes(5); % Argument of Perigee [rad]

P = 2*pi*sqrt(a^3/mu_E); % Keplerian period [s]

passLon = single(passLon);
% Discretise angular range in longitude
discLon = 3600;
Longitudes = single(linspace(0,360,discLon)');

% Matrix of longitudinal distance between pass and points on grid
normDeg = mod(bsxfun(@minus,Longitudes,repmat(passLon,length(Longitudes),1)),360);
normDegAbs = min(360-normDeg, normDeg);
% Find passes which pass within angular sensor width
inView = normDegAbs <= DlonA/2;

lat_max = min(inc_r, pi-inc_r)*(180/pi);

% Check for marginal pass at target latitude (Earth rotation + inclination)
%discTheta = 1000; % Discretisation

latRange = single(linspace(lat-theta,lat+theta)); % Latitude trace of sensor height
if any(latRange > lat_max)
    latRange(latRange > lat_max) = [];
    latRange = [latRange, lat_max];
end

if DscPass == 0
    nu = asin(sind(latRange)/sin(inc_r)) - AoP; % True Anomaly at Target Latitude [rad]
else
    latRange = fliplr(latRange);
    nu = pi - asin(sind(latRange)/sin(inc_r)) - AoP; % True Anomaly at Target Latitude [rad]
end
t_delta = (nu./(2*pi)) * P;

% Vary equivalent (elv or psi) sensor range with latitude
[theta_range,~,~,~] = calcTheta(latRange, coes, DscPass, f_p, f_e, flag, elv_psi);

% Compute equivalent (full) angular range (longitude) at equator
DlonA_range = 2*acosd((cosd(theta_range)-sind(latRange).^2)./cosd(latRange).^2) * (1+Flat);

%DlonA_range = DlonA_range * (0:
[~, midindx] = min(abs(latRange - lat));
t_delta = t_delta(1,:) - t_delta(1,midindx);

n = sqrt(mu_E/(a^3)); % Mean Motion [rad s^-1]
p = a*(1-ecc^2); % Semi-parameter [km]
% Compute drift in longitude crossing at equator [rad s^-1]
lonDrift = -omega_E -(3/2)*n*J2*(R_E/p)^2*cos(inc_r) + (3/32)*n*J2^2*(R_E^4/p^4)*cos(inc_r)*(12-4*ecc^2-(80+5*ecc^2)*sin(inc_r)^2);

% Longitude at Target Latitude [deg]
lonRange = atan2d(cos(AoP+nu)*sin(RAAN)+sin(AoP+nu)*cos(RAAN)*cos(inc_r), cos(AoP+nu)*cos(RAAN)-sin(AoP+nu)*sin(RAAN)*cos(inc_r))...
    + (nu/(2*pi))*P*lonDrift*(180/pi);
lonRangeAbs = lonRange - lonRange(midindx);
% if DscPass == 0
%     lonRangeAbs = lonRange - (lonRange(1) + (lonRange(end)-lonRange(1)/2)); % Absolute longitude trace of sensor width
% else
%     lonRangeAbs = (lonRange(1) + (lonRange(end)-lonRange(1)/2)) - lonRange; % Absolute longitude trace of sensor width
% end

% Vary distance of marginal view with target latitude, inclination, and FoV
% if (lat_max - lat) > DlonA
%     marginalDistance = DlonA;
% elseif (lat_max - lat) > DlonA/2
%     marginalDistance = DlonA*1.5;
% else
% end
% marginalDistance = (1/(lat_max - lat)) * 45;
marginalDistance = DlonA + max(min(abs(lonRangeAbs),360-abs(lonRangeAbs)));

% Define marginal view cases
marginalView = (normDegAbs <= marginalDistance);% - inView;
[idx(:,1),idx(:,2)] = find(marginalView == 1);
val = normDeg(marginalView == 1); % Distance between pass & grid points
lonMarg = mod(passLon(idx(:,2))',360); % Longitude of marginal passes
lonTrace = repmat(lonRangeAbs, length(idx), 1); % Array of dLon to test

% if DescendAscend == 1
%     nu2 = fliplr(pi-nu);
%     lonRange2 = atan2d(cos(AoP+nu2)*sin(RAAN)+sin(AoP+nu2)*cos(RAAN)*cos(inc_r), cos(AoP+nu2)*cos(RAAN)-sin(AoP+nu2)*sin(RAAN)*cos(inc_r));
%     lonRangeAbs2 = lonRange2 - lonRange2(discTheta/2); % Absolute longitude trace of sensor width
%     asc_dsc = mod(idx(:,2),2); % Asc = 1, Dsc = 1;
%     lonTrace = zeros(length(idx),discTheta);
%     lonTrace(asc_dsc == 1,:) = repmat(lonRangeAbs,sum(asc_dsc),1);
%     lonTrace(asc_dsc == 0,:) = repmat(lonRangeAbs2,length(idx)-sum(asc_dsc),1);
% end

lonMargMat = mod(bsxfun(@plus,lonTrace,lonMarg),360);

%lonDiff = mod(bsxfun(@minus, passLon(idx(:,2))',lonMargMat),360);
lonDiff = mod(bsxfun(@minus, mod(passLon(idx(:,2))'+val,360),lonMargMat),360);
lonDiff = min(360-lonDiff, lonDiff);

% Find if target is inside ellipse (x-h)^2/sMa^2 + (y-k)^2/sma^2 <= 1
a = bsxfun(@rdivide, lonDiff.^2,(0.5*DlonA_range).^2);
b = (lat-latRange).^2./theta_range.^2;
marginalPass = bsxfun(@plus, a, b);
marginalIdx = any(marginalPass <= 1 ,2);
t_delta = repmat(t_delta,length(marginalIdx),1);
t_view = t_delta;
t_view(marginalPass > 1) = 0;

indx = single(sub2ind(size(inView),idx(marginalIdx,1),idx(marginalIdx,2)));
inViewMarginal = false(size(inView));
inViewMarginal(indx) = 1;

t_view_min = single(zeros(size(inView)));
t_view_max = single(zeros(size(inView)));
t_view(t_view==0) = nan;
t_view_min(indx) = min(t_view(marginalIdx,:),[],2);
t_view_max(indx) = max(t_view(marginalIdx,:),[],2);

inView(indx) = 1;