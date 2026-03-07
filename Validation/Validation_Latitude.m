clc; clear;

% PURPOSE
% Computes revisit time for a fixed orbit while varying the
% target latitude. Useful to evaluate coverage performance
% across different regions of Earth.
%
% ORBIT PARAMETERS
% sma : semi-major axis [km above Earth radius]
% ecc : eccentricity
% inc : inclination [deg]
% RAAN: right ascension of ascending node [deg]
% AoP : argument of perigee [deg]
%
% SENSOR PARAMETERS
% elv : minimum elevation constraint [deg]
% psi : sensor half-cone angle (unused if elevation is defined)
%
% LATITUDE RANGE
% lat : vector of target latitudes [deg]
% Example:
% lat = 0:5:80
%
% CONSTELLATION SETTINGS
% Single satellite analysis:
% noSats   = 1
% noPlanes = 1
% relSpacing = 1
%
% FLAGS
% OneSideView     : 1 → sensor observes only one side of orbit
% DescendAscend   : 1 → accesses allowed on both ascending and descending passes
%
% TIME LIMITS
% dayLimit(1) : revisit search increment [days]
% dayLimit(2) : maximum analysis duration [days]
%
% OUTPUTS
% maxRevisit(lat) : maximum revisit time for each latitude [days]
%
% Converted result:
% maxRevisit_hr
%
% NOTES
% - Shows how revisit time varies from equator to high latitudes.
% - Particularly useful for analyzing polar or near-polar orbits.
% ================================================================

R_E = 6378.137;

sma = (R_E + 500);
ecc = 0;
inc = 97;
RAAN = 0;
AoP = 0;
coes = [sma*1e3, ecc, inc*(pi/180), RAAN, AoP];

dayLimit(1) = 1; % ListPasses Loop Control [days]
dayLimit(2) = 10; % Total Analysis Control [days]

OneSideView = 0; % Flag for one-side looking [0/1]
DescendAscend = 1; % Flag for access on both descending and ascending pass [0/1]

f_p = 0;
f_e = 1; elv = 30;
psi_elv = elv; 

lat = 0:5:80;

for i = 1:length(lat)
    tic
    [ maxRevisit(i), meanRevisit(i), ~] = RevisitCalc( coes, lat(i), f_p, f_e, 1, 1, 1, DescendAscend, OneSideView, dayLimit, psi_elv );
    toc
end

maxRevisit_hr = maxRevisit*24;

% psi = 0;
% [ MRT, ART ] = STK_Revisit_Fcn( sma, ecc, inc, elv, psi, lat);
% 
% MRT = MRT/3600;
% 
% Diff = MRT - maxRevisit;
% Error = Diff./MRT * 100;

% figure
% plot(

