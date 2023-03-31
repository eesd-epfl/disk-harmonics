
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
% Author of this file: Mahmoud S. Shaqfa and Gary Choi

%% Paths
% function Main_DH_roughness_projection
close all;
% clear; clc; tic
addpath('codes')
addpath('stlTools')
addpath('geodesic_domes')
addpath('input_geom/donor_meshes/')
addpath('analysis_results')
addpath('rec_output')
%% Reconstruction range for k \in [start, end]
starting_degree = 10 + 1;
ending_degree = 20 + 1;
scale_factor = 0.2;
touch_one_direction = true; % project vertices along one direction.
touched_direction = 3; % (1) -> x; (2) -> y; (3) -> z

plot_figures = true;
%% Load the analysis reults from HSH or SCH only
analysis_mat_file = 'rec_k_41_rocks_patch_rough.mat';

load(analysis_mat_file);
%% The donor mesh is the mesh that provides us with the shape and we project
% roughness on it extracted from the DH analyses.

donor_mesh = true; % if false it will load a hemisphere as a donor mesh
icosahedron_dome_refinement = 5; % Recommended refinement range 4-6
h_a = 1; % hemisphere dim along x
h_b = 1; % hemisphere dim along y

donor_mesh_name = "cylinder_ref2.stl";

if donor_mesh
    [v, f] = stlRead(donor_mesh_name);
else
%     [v, f] = icosphere(icosahedron_dome_refinement);
    [v, f] = solid_hemisphere(icosahedron_dome_refinement, h_a, h_b);
end
% v = [v(:,1) - mean(v(:,1)), v(:,2) - mean(v(:,2)), v(:,3) - mean(v(:,3))];
%% Specify (select) the group of vertices to be parametrized
% filter_query = strcat('selected_verts = find(v(:,3) < -0.1);'); % for a hemisphere
filter_query = strcat('selected_verts = find(v(:,3) > 0.5-0.1);'); % for a cylinder

eval(filter_query);
selected_groups = zeros(length(v(:, 1)), 1);
selected_groups(selected_verts) = 1;
fprintf('\nNumber of vertices chosen: %0.0f \n', length(selected_groups(selected_verts)))
%% Project roughness and render results
if plot_figures
    paraview_patch(v, f, selected_groups)
    title('Input donor mesh')
    view([50 20])
    colorbar off
end

if isempty(selected_verts)
    error('The selected vetices list must be specified.');
end

% Solve for \rho and \theta from the donor mesh patch
map = map_with_prescribed_disk(v, f, selected_verts);
if plot_figures
    plot_mesh(map, f)
    title('Input donor mesh parametrisation on a prescribed cap')
end
selected_map = map(selected_verts, :);
[phis, rhos] = cart2pol(selected_map(:,1), selected_map(:,2));
phis(phis<0) = phis(phis<0) + 2*pi(); % phis \in [0, 2*pi], rhos \in [0, 1]

k_range = [starting_degree, ending_degree];

partial_reconstruction = scale_factor .* real(DH_range_reconstruction_patch(k_range, qm_k, eigen_table, rhos, phis));

t = pi/4;
Rx = [1 0 0; 0 cos(t) -sin(t); 0 sin(t) cos(t)];
Ry = [cos(t) 0 sin(t); 0 1 0; -sin(t) 0 cos(t)];
Rz = [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1];

if touch_one_direction
    v(selected_verts, touched_direction) = v(selected_verts, touched_direction)...
        + partial_reconstruction(:, touched_direction, end);
else
    v(selected_verts, :) = v(selected_verts, :) + (Ry * partial_reconstruction(:, :, end)')';
end

if plot_figures
    paraview_patch(v, f, selected_groups)
    title('Output of the roughened mesh')
    colorbar off
end
%% Save the resulted STL file
stlWrite(strcat('roughness_output/rec_l_', num2str(starting_degree), '_to_',...
        num2str(ending_degree), '_DH_', num2str(h_a), '_', num2str(h_b), surface_name)...
        , f, v,'mode','ascii')
save(strcat('analysis_results/roughness_rec_k_', num2str(max_degree), '_', surface_name_wo_format,...
    '.mat'))
% end