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
% Author of this file: Gary Choi

function map = DEM_disk(v,f,population,initial_map,epsilon,step_ratio)
% main program for disk/elliptic density-equalizing map 

% Input:
% v: m x 2 or m x 3 vertex coordinates
% f: k x 3 triangulation
% population: k x 1 positive array
% initial_map: m x 2 coordinates
% (Optional) epsilon: stopping parameter (default = 0.001)
% (Optional) step_ratio: rescale the step size (default = 1)

% output:
% map: m x 2 density-equalizing flattening map 
% Remark: the area of the result will be the same to input mesh


% default stopping parameter
if ~exist('epsilon','var')
    epsilon = 1e-3;
end

% rescale the step size by a multiple
if ~exist('step_ratio','var')
    step_ratio = 1;
end

% find boundary
TR = TriRep(f,v); bdy_e = freeBoundary(TR); bdy_v = bdy_e(:,1); 
bdy_f = zeros(length(bdy_e),1);
for i = 1:length(bdy_e)
    bdy_f(i) = find(sum((f == bdy_e(i,1)) + (f==bdy_e(i,2)),2) == 2);
end

area_original = sum(face_area(f,v));

r = initial_map(:,1:2);
r = r * sqrt(area_original / sum(face_area(f,r)));

r_bdy = r(bdy_v,:);

% density
rho_f = population ./ face_area(f,r); 
rho_v = f2v_area(r,f)*rho_f; 

step = 1;

dt = min([min(rho_v)/mean(rho_v),mean(rho_v)/max(rho_v)])*area_original*step_ratio;

disp('Step     std(rho)/mean(rho)');  
disp([num2str(step), '        ',num2str(std(rho_v)/mean(rho_v))]);
    
rho_v_error = zeros(1,500);
rho_v_error(step) = std(rho_v)/mean(rho_v);

while rho_v_error(step) > epsilon
    % update rho 
    L = laplace_beltrami(r,f);
    A = lumped_mass_matrix(r,f);
    rho_v_temp = (A - dt*L) \ (A*rho_v);
    
    if sum(sum(isnan(rho_v_temp)))~=0
        break;
    end

    % update gradient 
    grad_rho_temp_f = compute_gradient(r,f,rho_v_temp);
    grad_rho_temp_v = f2v_area(r,f)*grad_rho_temp_f;
    
    % projection
    bdy_normal1 = [-(r(bdy_e(:,1),2) - r(bdy_e(:,2),2)),  (r(bdy_e(:,1),1) - r(bdy_e(:,2),1))];
    bdy_normal1 = [bdy_normal1(:,1) ./ sqrt(bdy_normal1(:,1).^2 + bdy_normal1(:,2).^2) , bdy_normal1(:,2) ./ sqrt(bdy_normal1(:,1).^2 + bdy_normal1(:,2).^2)];
    bdy_normal2 = [-(r(bdy_e([end,1:end-1],1),2) - r(bdy_e([end,1:end-1],2),2)),  (r(bdy_e([end,1:end-1],1),1) - r(bdy_e([end,1:end-1],2),1))];
    bdy_normal2 = [bdy_normal2(:,1) ./ sqrt(bdy_normal2(:,1).^2 + bdy_normal2(:,2).^2) , bdy_normal2(:,2) ./ sqrt(bdy_normal2(:,1).^2 + bdy_normal2(:,2).^2)];
    bdy_normal = ( bdy_normal1 + bdy_normal2 )/2;
    
    grad_rho_temp_v(bdy_v,:) = grad_rho_temp_v(bdy_v,:) - ...
        [sum(bdy_normal.*grad_rho_temp_v(bdy_v,:),2).*bdy_normal(:,1), ...
        sum(bdy_normal.*grad_rho_temp_v(bdy_v,:),2).*bdy_normal(:,2)];
    
    % update displacement 
    dr = [grad_rho_temp_v(:,1) ./ rho_v_temp, grad_rho_temp_v(:,2) ./ rho_v_temp];
    r_temp = r - dr*dt;
    r = r_temp; 
    
    % update other parameters 
    rho_v = rho_v_temp;
    step = step + 1;

    rho_v_error(step) = std(rho_v)/mean(rho_v);
    disp([num2str(step), '        ',num2str(std(rho_v)/mean(rho_v))]);

    if step > 500
        break;
    end
end

% a projection to correct any boundary shape mismatch due to numerical approximation error
proj_bdy = distance2curve(r_bdy,r(bdy_v,:));
r(bdy_v,:) = proj_bdy;

% finally normalize the area
map = r * sqrt(area_original / sum(face_area(f,r)));

fprintf('Completed.\n');