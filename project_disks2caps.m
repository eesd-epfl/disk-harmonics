%*************************************************************************%
% All rights reserved (C) to the authors: Mahmoud SHAQFA, Gary CHOI,      %
% Guillaume ANCIAUX & and Katrin BEYER                                    %
%                                                                         %
% M. Shaqfa Contact:                                                      %
% Department of Mechanical Engineering, Massachusetts Institute of        %
% Technology (MIT)                                                        %
% Cambridge, MA, USA                                                      %
%               Email: mshaqfa@mit.edu                                    %
%                                                                         %
% G. Choi Contact:                                                        %
% Department of Mathematics, Massachusetts Institute of Technology (MIT)  %
% Cambridge, MA, USA                                                      %
%               Email: ptchoi@mit.edu                                     %
%                                                                         %
% G. Anciaux Contact:                                                     %
%               Email: guillaume.anciaux@epfl.ch                          %
%                                                                         %
% K. Beyer Contact:                                                       %
%               Email: katrin.beyer@epfl.ch                               %
%*************************************************************************%
% This code includes implementations for:                                 %
%				- Disk harmonics expansion for parametric surfaces        %
% This code is part of the paper: "Disk Harmonics for Analysing Curved    %
% and Flat Self-affine Rough Surfaces and the Topological                 %
% Reconstruction of Open Surfaces"                                        %
%                                                                         %
%*************************************************************************%
% This library is free software; you can redistribute it and/or modify	  %
% it under the terms of the GNU Lesser General Public License as published%
% by the Free Software Foundation; either version 2.1 of the License, or  %
% (at your option) any later version.                                     %
%                                                                         %
% This library is distributed in the hope that it will be useful,         %
% but WITHOUT ANY WARRANTY; without even the implied warranty of          %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    %
% See the GNU Lesser General Public License for more details.        	  %
% You should have received a copy of the GNU Lesser General Public License%
% along with this library; if not, write to the Free Software Foundation, %
% Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA       %
%*************************************************************************%
% Author of this file: Mahmoud S. Shaqfa, Gary Choi

function map = project_disks2caps(v, r, theta_c)

area_sphericalcap = 2*pi*r*(1-cos(theta_c));

area_original_disk = pi * (max(sqrt(v(:,1).^2 + v(:,2).^2)))^2;

scale = sqrt(area_sphericalcap/area_original_disk); % should be with sqrt??

% Normalize the vertices for a unit disk
v = [v(:,1) - mean(v(:,1)), v(:,2) - mean(v(:,2)), v(:,3) - mean(v(:,3))];
v(:, 1:2) = v(:, 1:2) ./ (max(sqrt(v(:,1).^2 + v(:,2).^2)));

p = v;
p = p*sqrt(2*(1-cos(theta_c))); % Scale it for a spherical cap
% Lambert equi-area projection
v_rec = -[sqrt(1-(p(:,1).^2+p(:,2).^2)/4).*p(:,1), ...
    sqrt(1-(p(:,1).^2+p(:,2).^2)/4).*p(:,2), ...
    -1+(p(:,1).^2+p(:,2).^2)/2];
x = v_rec(:, 1); y = v_rec(:, 2); z = v_rec(:, 3);

%%
% Scale the sphere if needed
x = x.*r; y = y.*r; z = z.*r;

% add the noise along the normal direction
noise = v(:,3) .*scale;
x = x.*(1+noise);
y = y.*(1+noise);
z = z.*(1+noise);
map = [x,y,z];

end