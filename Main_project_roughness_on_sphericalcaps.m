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
% Author of this file: Mahmoud S. Shaqfa and Gary Choi

close all;
clear; clc; tic
addpath('codes')
addpath('stlTools')
addpath('input_geom')
addpath('distmesh')


surface_name = "disk_rec_artificial_H_0.95.stl";

[v,f] = stlRead(surface_name);

theta_c = 10;
theta_c = deg2rad(theta_c);

r=3;
area_sphericalcap = 2*pi*r*(1-cos(theta_c));
area_original_disk = pi * (max(sqrt(v(:,1).^2 + v(:,2).^2)))^2;

scale = sqrt(area_sphericalcap/area_original_disk);


% Normalize the vertices for a unit disk
v = [v(:,1) - mean(v(:,1)), v(:,2) - mean(v(:,2)), v(:,3) - mean(v(:,3))];
v(:, 1:2) = v(:, 1:2) ./ (max(sqrt(v(:,1).^2 + v(:,2).^2)));

paraview_patch(v, f)

convert = true;
projection_type = 2;

if convert
    % use stereographic projection here
    if projection_type == 1
        % can also use the Lambert projection or other method
        X = v(:, 1); Y = v(:, 2);
        X_normalized = X;
        Y_normalized = Y;
        x = 2*X_normalized./(1+X_normalized.^2+Y_normalized.^2);
        y = 2*Y_normalized./(1+X_normalized.^2+Y_normalized.^2);
        z = (1-X_normalized.^2-Y_normalized.^2)./(1+X_normalized.^2+Y_normalized.^2);
    else
        p = v;
        p = p*sqrt(2*(1-cos(theta_c))); % Scale it for a spherical cap
        % Lambert equi-area projection
        v_rec = -[sqrt(1-(p(:,1).^2+p(:,2).^2)/4).*p(:,1), ...
            sqrt(1-(p(:,1).^2+p(:,2).^2)/4).*p(:,2), ...
            -1+(p(:,1).^2+p(:,2).^2)/2];
        x = v_rec(:, 1); y = v_rec(:, 2); z = v_rec(:, 3);
    end
    
    % Scale the sphere if needed
    x = x.*r; y = y.*r; z = z.*r;
    % add the noise along the normal direction
    noise = v(:,3) .*scale;
    x = x.*(1+noise);
    y = y.*(1+noise);
    z = z.*(1+noise);
    rec_v = [x,y,z];
    
    figure;
    scatter3(x,y,z);
    axis equal tight
    xlabel('x');
    ylabel('y');
    colorbar
    title('Noise')
    stlWrite(strcat('input_geom/', 'sphericalcap_th_', num2str(rad2deg(theta_c)), '_r_', num2str(r),'_',...
        surface_name), f, rec_v,'mode','ascii')
else
    disp("Not convereted")
end
