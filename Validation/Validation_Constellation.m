clear
clc

% PURPOSE
% Computes revisit time for different satellite constellation
% configurations. Each index i of the different vectors represents a different constellation.
%
% PARAMETERS TO MODIFY
% sma = 6378.137 + [700 1100 1500];  : orbital altitude vector [km above Earth radius]
% ecc = zeros(1,10);                 : eccentricity
% inc = [90 86 96]                   : inclination [deg]
% elv = [0 10 20];                   : minimum elevation angle constraint [deg]
% lat = zeros(1,3);                  : target latitude [deg]
% psi = [];                          : sensor half cone angle [deg] - not used here
%
% noSats = [3 3 3];                  : total number of satellites in constellation
% noPlanes = [3 3 3];                : number of orbital planes
% relSpacing  = [0 0 1];             : relative satellite spacing between planes
%
% SENSOR PARAMETERS
% elv        : elevation mask angle (ground station constraint)
% psi        : sensor half-cone angle (if used instead of elevation)
%
% FLAGS
% OneSideView     : 1 → sensor only looks on one side of orbit
% DescendAscend   : 1 → accesses allowed on ascending AND descending passes
%
% TIME LIMITS
% dayLimit(1) : revisit search step [days]
% dayLimit(2) : maximum analysis duration [days]
%
% OUTPUTS
% maxRevisit  : maximum revisit time [days]
% meanRevisit : average revisit time [days]
%
% Converted results:
% maxRevisit_hr
% maxRevisit_min
% maxRevisit_sec
%
% NOTES
% - Each index i corresponds to a different constellation configuration.
% - Useful for comparing Walker-like constellations with different
%   altitudes, inclinations, and satellite distributions.
% ================================================================


sma = 6378.137 + [700 1100 1500];
ecc = zeros(1,10);
inc = [90 86 96];
elv = [0 10 20];
lat = zeros(1,3);
psi = [];

noSats = [3 3 3]; 
noPlanes = [3 3 3]; 
relSpacing = [0 0 1];

OneSideView = 0; % Flag for one-side looking [0/1]
DescendAscend = 1; % Flag for access on both descending and ascending pass [0/1]

dayLimit(1) = 1; % ListPasses Loop Control [days]
dayLimit(2) = 5; % Total Analysis Control [days]

f_p=~isempty(psi); % Psi flag
f_e=~isempty(elv);  % Elv flag

% f_p=exist('psi','var'); % Psi flag
% f_e=exist('elv','var'); % Elv flag
psi_elv = elv; 

for i = 1:length(sma)
    coes = [sma(i)*1e3, ecc(i) ,inc(i)*(pi/180), 0, 0];
    tic
    [ maxRevisit(i,1), meanRevisit(i,:) ] = RevisitCalc( coes, lat(i), f_p, f_e, noSats(i), noPlanes(i), relSpacing(i), DescendAscend, OneSideView, dayLimit, psi_elv(i) );
    toc
end

maxRevisit_hr = maxRevisit*24;
maxRevisit_min = maxRevisit*24*60;
maxRevisit_sec = maxRevisit*24*60*60;
