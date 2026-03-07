clear
clc

% PURPOSE
% Evaluates revisit performance for multiple independent orbit cases
% using a single satellite.
%
% PARAMETERS TO MODIFY
% sma : orbital altitude [km above Earth radius]
% ecc : eccentricity
% inc : inclination [deg]
% elv : minimum elevation constraint [deg]
% lat : target latitude [deg]
%
% Each index i represents a different orbit configuration.
%
% CONSTELLATION SETTINGS
% noSats   = 1
% noPlanes = 1
% relSpacing = 0
%
% SENSOR PARAMETERS
% elv : minimum elevation angle constraint [deg]
% psi : sensor half-cone angle (leave empty if elevation used)
%
% FLAGS
% OneSideView     : 1 → sensor only observes on one side
% DescendAscend   : 1 → accesses allowed on both ascending and descending passes
% SSO             : Sun-synchronous orbit flag (not used here)
%
% TIME LIMITS
% dayLimit(1) : revisit search increment [days]
% dayLimit(2) : maximum analysis duration [days]
%
% OUTPUTS
% maxRevisit  : maximum revisit time [days]
% meanRevisit : mean revisit time [days]
%
% Converted results:
% maxRevisit_hr
% maxRevisit_min
% maxRevisit_sec
%
% NOTES
% - Useful for orbit design trade studies.
% - Allows comparison of altitude, inclination, and sensor elevation
%   constraints for a single spacecraft.
% ================================================================

sma = 6378.137 + [400 400 400 400 800 800 800 800 550 700];
ecc = zeros(1,10);
inc = [20 20 60 60 20 20 60 60 97.59 98.19];
elv = [10 40 10 40 10 40 10 40 20 30];
psi = [];
lat = zeros(1,10);

noSats = 1; noPlanes = 1; relSpacing = 0;

OneSideView = 0; % Flag for one-side looking [0/1]
DescendAscend = 1; % Flag for access on both descending and ascending pass [0/1]
SSO = 0; % SSO orbit flag
dayLimit(1) = 1; % ListPasses Loop Control [days]
dayLimit(2) = 5; % Total Analysis Control [days]

f_p=~isempty(psi); % Psi flag
f_e=~isempty(elv);  % Elv flag

psi_elv = elv; 

for i = 1:length(sma)
    coes = [sma(i)*1e3, ecc(i) ,inc(i)*(pi/180), 0, 0];
    tic
    [ maxRevisit(i,1), meanRevisit(i,:), flag ] = RevisitCalc( coes, lat(i), f_p, f_e, noSats, noPlanes, relSpacing, DescendAscend, OneSideView, dayLimit, psi_elv(i) );
    toc
end

maxRevisit_hr = maxRevisit*24;
maxRevisit_min = maxRevisit*24*60;
maxRevisit_sec = maxRevisit*24*60*60;
