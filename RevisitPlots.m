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
    % Contour / Line plot of Maximum Revisit Time [day] vs altitude [km] and inclination [deg]
    if SSO == 0
        % Prepare Z for indexing: ensure MaxRevisit dims are [nInc, nHa, nPsi] or similar
        Z = squeeze(MaxRevisit(:,:,1)); % dimensions may be [nInc x nHa] or [1 x nHa] or [nInc x 1]

        nInc = numel(inc);
        nHa  = numel(h_a);

        % Both vary -> contour
        if nInc >= 2 && nHa >= 2
            figure(4)
            % contour expects Z to be size length(Y) x length(X) where X is h_a, Y is inc
            % If Z currently is [nHa x nInc], transpose it
            if isequal(size(Z), [nHa, nInc])
                Z = Z';
            end
            contour(h_a/1e3, inc, Z, 'fill', 'on')
            xlabel('Altitude [km]')
            fnt = ylabel('Inclination [deg]');
            c = colorbar;
            c.Label.String = 'Maximum Revisit Time [day]';
            c.Label.FontSize = get(fnt, 'fontsize');
            if flag == 1
                title({['Contour plot of Maximum Revisit Time for'];[' varying orbit altitude and inclination'];['\epsilon = ', num2str(elv), ' deg, \phi = ', num2str(lat), ' deg']})
            elseif flag == 0 && numel(psi) == 1
                title({['Contour plot of Maximum Revisit Time for'];[' varying orbit altitude and inclination'];['\psi = ', num2str(psi), ' deg, \phi = ', num2str(lat), ' deg']})
            end

            % Only altitude varies -> 1-D plot vs altitude
        elseif nHa >= 2 && nInc <= 1
            figure(4)
            plot(h_a/1e3, reshape(Z, 1, []))
            xlabel('Altitude [km]'), ylabel('Maximum Revisit Time [day]')

            % Only inclination varies -> 1-D plot vs inclination
        elseif nInc >= 2 && nHa <= 1
            figure(4)
            plot(inc, reshape(Z, [], 1))
            xlabel('Inclination [deg]'), ylabel('Maximum Revisit Time [day]')

        else
            warning('Not enough points to plot altitude/inclination surface. Need at least two values in one dimension.');
        end
    end
    figure (5)
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
