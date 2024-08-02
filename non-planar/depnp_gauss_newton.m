function [x_opt,err,depth,depth_Xc]=depnp_gauss_newton(CB,Alpha,impts,Alphr,imptPart,x0)

% GAUSS_NEWTON  
%
%       K 12x4:  Matrix of kernel vectors
%       Cw 4x3:  Coordinates of the control points in world reference
%       Beta0: initial estimation of Beta's
%
%
% Copyright (C) <2007>  <Francesc Moreno-Noguer, Vincent Lepetit, Pascal Fua>
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the version 3 of the GNU General Public License
% as published by the Free Software Foundation.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.       
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%
% Francesc Moreno-Noguer, CVLab-EPFL, October 2007.
% fmorenoguer@gmail.com, http://cvlab.epfl.ch/~fmoreno/ 

current_x=x0;
errStop = 1e-6;
n_iterations=10; %max number of iterations. Usually 4-5 is enought.
error = zeros(n_iterations,1);
for k=1:n_iterations

    [J,e,depth,depth_Xc]=depnp_GN_Ab(CB,Alpha,impts,Alphr,imptPart,current_x);  
    dbeta=(J'*J)\J'*e;
    current_x=current_x-dbeta;
    error(k)=e'*e;
    if ((k>2)&&(abs(error(k)-error(k-1))<errStop) )
        break;
    end
end

x_opt=current_x;
err=error(end);

end


