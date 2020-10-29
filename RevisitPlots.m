% RevisitPlots provides examples for plotting the revisit time metrics
% associated with the RevisitGrid script.
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

% XY Plot of Maximum Revisit Time [day] vs. altitude
figure (1)
plot(h_a/1e3, squeeze(MaxRevisit(1,:,1)))
xlabel('Altitude [km]'), ylabel('Maximum Revisit Time [day]')
ylim([0, min([max(max((MaxRevisit))),dayLimit(2)])]); 

% Contour plot of Maximum Revisit Time [day] for varying psi angle [deg] and altitude [km]
if flag == 0 && numel(psi)>=2 
    figure (2)
    contour(h_a/1e3, psi, squeeze(MaxRevisit(1,:,:))', 'fill', 'on')
    xlabel('Altitude [km]'), fnt = ylabel('\psi [deg]');
    c = colorbar;
    c.Label.String = 'Maximum Revisit Time [day]';
    c.Label.FontSize = get(fnt, 'fontsize');
    title({['Contour plot of Maximum Revisit Time for'];[' varying orbit altitude and Half-Cone angle']})
end

% Contour plot of Maximum Revisit time [day] according to altitude [km] and inclination [deg]
if SSO == 0
    figure (3)
    plot(inc, squeeze(MaxRevisit(:,1,1)))
    xlabel('Inclination [km]'), ylabel('Maximum Revisit Time [day]')
    figure (4)
    contour(h_a/1e3, inc, squeeze(MaxRevisit(:,:,1)), 'fill', 'on')
    xlabel('Altitude [km]'), fnt = ylabel('Inclination [deg]') ;
    c = colorbar;
    c.Label.String = 'Maximum Revisit Time [day]';
    c.Label.FontSize = get(fnt, 'fontsize');
    if flag == 1
        title({['Contour plot of Maximum Revisit Time for'];[' varying orbit altitude and inclination'];['\epsilon = ', num2str(elv), ' deg, \phi = ', num2str(lat), ' deg']})
    elseif flag == 0 && numel(psi) == 1
        title({['Contour plot of Maximum Revisit Time for'];[' varying orbit altitude and inclination'];['\psi = ', num2str(psi), ' deg, \phi = ', num2str(lat), ' deg']})
    end
end

% Contour plot of Maximum Revisit time [day] according to altitude [km] and inclination [deg], angle of view (psi) [deg] is fixed
if SSO == 0 && flag == 0 && numel(psi) == 1 
    figure (1)
    contour(h_a/1e3, inc, squeeze(MaxRevisit(:,:,1)), 'fill', 'on')
    xlabel('Altitude [km]'), ylabel('Inclination [deg]')
    c = colorbar;
    c.Label.String = '\fontsize{10.5} Maximum Revisit Time [day]';
title({['Contour plot of Maximum Revisit Time for'];[' varying orbit altitude and inclination'];['\psi = ', num2str(psi), ' deg, \phi = ', num2str(lat), ' deg']})
end
