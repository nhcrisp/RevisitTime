function [MaxRevisit] = ReShape(maxRevisit, flag, SSO, inclination, altitude, psi)
%RESHAPE Reshaping Max Revisit vector into a n-dimensions matrix
%
% Inputs:
%       maxRevisit      Vector of maximum revisit time [days]
%       flag            Flag for Elv override
%       SSO             Flag for SSO
%       inclination     Vector of inclinations
%       altitude        Vector of altitudes
%       psi             Vector of FoV
%
% Outputs:
%       MaxRevisit      Matrix of maximum revisit time [days]
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
j = 1;

if flag == 0 && SSO == 0
    for ind_inc = 1:inclination
        for ind_alt = 1:altitude
            for ind_psi = 1:psi
                MaxRevisit(ind_inc, ind_alt, ind_psi) = maxRevisit(j);
                j=j+1;
            end
        end
    end
elseif flag == 1 && SSO == 0
    for ind_inc = 1:inclination
        for ind_alt = 1:altitude
            MaxRevisit(ind_inc, ind_alt, 1) = maxRevisit(j);
            j=j+1;
        end
    end
elseif flag == 0 && SSO == 1
    for ind_alt = 1:altitude
        for ind_psi = 1:psi
            MaxRevisit(1,ind_alt, ind_psi) = maxRevisit(j);
                j=j+1;
        end
    end
elseif flag == 1 && SSO == 1
    for ind_alt = 1:altitude
        MaxRevisit(1,ind_alt) = maxRevisit(j);
        j=j+1;
    end
end

end
    
