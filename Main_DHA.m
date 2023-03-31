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
% Author of this file: Mahmoud S. Shaqfa

% % This is the main script interface for the Disk Harmonic Analysis
%% Paths
close all;
clear; clc; tic
addpath('codes')
addpath('stlTools')
addpath('input_geom')
addpath('distmesh')
%% Inputs
fprintf('The analysis started at: %s\n', datestr(now,'HH:MM:SS.FFF'))

% To appply spherical cap parametrization step (true)
parametrization = true;
python_solver = false;

surface_name = "Matterhorn_mode.stl";

% mapping_type = "conformal"; % This option is not very good for reconstruction.
mapping_type = "area_preserving";


% Note: these roots have been imported from Python precomputed tables of
% special functions. They're more reliable and accurate than the computed
% ones here.
para_surface_path = ''; % if the parametrization already computed
solve_roots = false; % To solve for the roots if not pre-calculated!
eigen_table_path = "eigenvalues_k_200.mat";

max_degree = 5 + 1; % +1 to account for the 0-degree.
max_reconstruction_degree = 70+ 1;
truncation_degree_fit = max_reconstruction_degree; % For fitting the fractal dimension
truncation_degree_fit = min(truncation_degree_fit, max_degree);
max_reconstruction_degree = min(max_reconstruction_degree, max_degree);

solver_increment_diag = 0.1;
solver_increment_tri = 0.1;
iterations_limit = 100;

% Math grid reconstruction
reconstruction_resolution = 100;
math_grid_reconstruction = false;

% Geodesic dome STL reconstruction
recursive_reconstruction = true;
recursive_increments = 10;

% option 1 gives you reconstruction by edge length and 2 gives you
% reconstruction by number of even elements used.
reconstruction_mesh_option = 1;
% The desired length of each edge on a unit disk from (0.035 < edge_length < 1.0)
edge_length = 0.0150;
% The resolution or the number of vertices must be an even number!
resolution = 50*2;

plot_figures = true;
%% Spherical Cap Parametrization

% Load the STL file
[v,f] = stlRead(surface_name);

% Normalize the location of the model (optional to uncomment)
% v = [v(:,1) - mean(v(:,1)), v(:,2) - mean(v(:,2)), v(:,3) - mean(v(:,3))];
% v(:,3) = v(:,3) + ones(length(v(:,3)), 1).*5;

if parametrization
    % ensure there is no triangle with all three vertices lying on the boundary
    % (such a triangle will be highly squeezed under the parameterization as
    % all three vertices will be constrained)
    v_old = v;
    f_old = f;
    [v,f] = clean_mesh(v_old,f_old,1);
    while length(v) ~= length(v_old)
        v_old = v;
        f_old = f;
        [v,f] = clean_mesh(v_old,f_old,1);
    end
    if plot_figures
        height_map = v(:, 3);
        paraview_patch(v, f, height_map)
        view([50 20])
        title('Input surface')
    end
    if mapping_type == "conformal"
        map = disk_conformal_map(v,f);
    elseif mapping_type == "area_preserving"
        map = disk_density_equalizing_map(v,f);
    end
    map(:, 3) = 0; % Add the z-column with zeros
    if plot_figures
        paraview_patch(map, f, height_map)
        view([0 90])
        title('Planar disk parameterization')
        axis equal tight on
    end
    % Save parametrised mesh
    try
        stlWrite(strcat('rec_output/para_output_', surface_name), f, map,'mode','ascii')
    catch
        stlWrite(strjoin(['rec_output/para_output_', surface_name], ""), f, map,'mode','ascii')
    end
else
    [map,f] = stlRead(para_surface_path);
end
fprintf("\n\n Verts number: %3.0f", length(v(:, 1)))
%% Finding the Sturm–Liouville eigenvalues for the spherical cap basis.
% This step should be ideally run from the "Sturm_Liouville_eigenvalues.m"
% file and not here so it can be adjusted and calibrated.
if solve_roots
    eigen_table = Sturm_Liouville_eigenvalues(max_degree, "neumann",...
        solver_increment_diag, solver_increment_tri, iterations_limit);
else
    % Load the table
    eigen_table = load(eigen_table_path);
    eigen_table = eigen_table.eigen_table;
end
%% The analysis step (DHA)
[phis, rhos] = cart2pol(map(:,1), map(:,2));
phis(phis<0) = phis(phis<0) + 2*pi(); % phis \in [0, 2*pi], rhos \in [0, 1]

D_mat = disk_harmonic_basis(max_degree, eigen_table, rhos, phis); % Disk harmonics basis

if python_solver
    fprintf("\n\nPython solver initiated...\n")
    save("mat2py_vars.mat", "D_mat", "v", '-v7.3')
    systemCommand = ['python3.8 python_solver.py']; % Change the version 3.8 if needed!
    [status, results] = system(systemCommand);
    if status == 0
        disp(results)
        load("py2mat_vars.mat");
        qm_k = qm_k{1,1};
        delete("mat2py_vars.mat")
        delete("py2mat_vars.mat")
    else
        disp("Python solver failed!!")
    end
else
    t1 = toc;
    qm_k = D_mat\v;
    t2 = toc-t1;
    fprintf("\n\nMatlab's Pseudo solver time: %3.6f seconds\n", t2)
end
%% Shape descriptors and the fractal dimension
% Shape descriptors (2-norm) for frequency accumulates at a certain
% frequency degree..
Dl = zeros([3, max_degree+1]);
for  k = 1:3
    for l = 1:max_degree
        for m = -l:1:l
            Dl(k, l) = Dl(k, l) + (real(qm_k(l^2 + l + m + 1, k)))^2 + (imag(qm_k(l^2 + l + m + 1, k)))^2;
        end
        Dl(k, l) = sqrt(Dl(k, l));
    end
    Dl(k, :) = Dl(k, :)/Dl(k, 1);
end
D0 = zeros([4, max_degree+1]);
for  k = 1:3
    m = 0;
    for l = 1:max_degree
        D0(k, l) = D0(k, l) + (real(qm_k(l^2 + l + m + 1, k)))^2 + (imag(qm_k(l^2 + l + m + 1, k)))^2;
        D0(k, l) = sqrt(D0(k, l));
        if k == 3
            D0(4, l) = eigen_table(l+1, abs(m)+1);
        end
    end
    D0(k, :) = D0(k, :)/D0(k, 1);
end

if plot_figures
    figure
    subplot(3,3,1)
    x_temp = 1:length(Dl(1, :));
    loglog(x_temp(2:end), Dl(1, 2:end), 'LineWidth', 2)
    title('The shape descriptors')
    xlabel('Freqency index (k)')
    ylabel('Normalized amplitude (Dx)')
    grid on
    
    subplot(3,3,2)
    loglog(x_temp(2:end), Dl(2, 2:end), 'LineWidth', 2)
    title('The shape descriptors')
    xlabel('Freqency index (k)')
    ylabel('Normalized amplitude (Dy)')
    grid on
    
    subplot(3,3,3)
    loglog(x_temp(2:end), Dl(3, 2:end), 'LineWidth', 2)
    title('The shape descriptors')
    xlabel('Freqency index (k)')
    ylabel('Normalized amplitude (Dz)')
    grid on
    
    subplot(3,3,4:6)
    loglog(x_temp(2:end), sqrt(Dl(1, 2:end).^2 + ...
        Dl(2, 2:end).^2 + Dl(3, 2:end).^2), 'LineWidth', 2)
    title(['The shape descriptors'])
    xlabel('Freqency index (k)')
    ylabel('Normalized amplitude (Dr)')
    grid on
    
    subplot(3,3,7:9)
    loglog(x_temp(2:truncation_degree_fit), (Dl(1, 2:truncation_degree_fit).^2 + ...
        Dl(2, 2:truncation_degree_fit).^2 + Dl(3, 2:truncation_degree_fit).^2), 'LineWidth', 2)
    title(['Surface Energy'])
    xlabel('Freqency index (k)')
    ylabel('Energy of amplitude (Dr^{2})')
    grid on
    hold on
    loglog(x_temp(2:truncation_degree_fit), (D0(1, 2:truncation_degree_fit).^2 + ...
        D0(2, 2:truncation_degree_fit).^2 + D0(3, 2:truncation_degree_fit).^2), 'LineWidth', 2)
    hold off
end
surface_name_wo_format = split(surface_name,'.');
surface_name_wo_format = surface_name_wo_format{1};
save(['rec_output/analysis_k_', num2str(max_degree), ...
    '_', surface_name_wo_format, '.mat'], 'Dl')

%% Estimate the wavelengths by FDEC (This need validation not correct so far for DHA)
FDEC = qm_k(2:4, :);

l11 = eigen_table(1+1, abs(1)+1);
l10 = eigen_table(1+1, abs(0)+1);
N_11 = 0.5 * (1 - (abs(1)/l11)^2) * besselj(abs(1), l11)^2;
N_11 = 1/sqrt(N_11);
N_11 = N_11 * (1/sqrt(2*pi));
N_10 = 0.5 * (1 - (abs(0)/l10)^2) * besselj(abs(0), l10)^2;
N_10 = 1/sqrt(N_10);
N_10 = N_10 * (1/sqrt(2*pi));

A = [N_11 .* (FDEC(:, 1) - FDEC(:, 3)), N_11 .* 1i .* (FDEC(:, 1) + FDEC(:, 3)), N_10 .* FDEC(:, 2)];

[vv,dd] = eig(A*A');
d = sqrt(diag(abs(dd)));
omega = (2*pi/k) * sqrt( ((0.5*d(2))^2 + (0.5*d(3))^2)/2);
fprintf("\n\n The dimensions of FDEC: a = %3.3f, b = %3.3f, c = %3.3f.\n",...
    d(3), d(2), d(1))
x_omega = 1:(max_degree-2);
y_omega = (2*pi) * sqrt(((0.5*d(2))^2 + (0.5*d(3))^2)/2) ./ x_omega;
y_omega_2 = 2*sqrt((0.5*d(1))^2 + (0.5*d(2))^2 + (0.5*d(3))^2) ./ (x_omega - 0.25);

if plot_figures
    figure
    loglog(x_omega, y_omega, 'LineWidth', 2)
    xlabel('Freqency index (k)')
    ylabel('Wavelength (\omega_{k})')
    title('The corresponding wavelengths')
end
fprintf("\n\n The largest angular wavelength = %3.4f, and the smallest = %3.4f.\n",...
    max(y_omega), min(y_omega))
fprintf("\n\n The largest radial wavelength = %3.4f, and the smallest = %3.4f.\n",...
    max(y_omega_2), min(y_omega_2))

%% The reconstruction step
fprintf("\n\n The reconstruction part has started...")
if math_grid_reconstruction
    reconstructed = real(disk_harmonic_math_reconstruction(max_reconstruction_degree, ...
        eigen_table, reconstruction_resolution, qm_k));
    figure;
    surf(reconstructed(:, :, 1), reconstructed(:, :, 2), reconstructed(:, :, 3));
    axis equal tight off
    title(['Reconstructed surface', ' with K = ', num2str(max_reconstruction_degree)]);
    view([50 20])
    % Define colormap
    red_color = zeros(1, 255); green_color = zeros(1, 255); blue_color = zeros(1, 255);
    red_color(1:127) = linspace(42,220,127); red_color(128:255) = linspace(220,174,127+1);
    green_color(1:127) = linspace(63,220,127); green_color(128:255) = linspace(220,0,127+1);
    blue_color(1:127) = linspace(181,220,127); blue_color(128:255) = linspace(220,22,127+1);
    ParaviewMap = [red_color', green_color', blue_color']./255;
    colormap(ParaviewMap);
    s.CData = reconstructed(:, :, 3);
    colorbar
else
    [v_rec, f_rec] = uniform_disk_grid(edge_length, resolution,...
        reconstruction_mesh_option);
    v_rec(:, 3) = 0;
    [rec_phis, rec_rhos] = cart2pol(v_rec(:,1), v_rec(:,2));
    rec_phis(rec_phis<0) = rec_phis(rec_phis<0) + 2*pi(); % phis \in [0, 2*pi]
    rec_basis_func = disk_harmonic_basis(max_degree, eigen_table, rec_rhos, rec_phis);
    v_rec = real(rec_basis_func * qm_k);
    try
        stlWrite(strcat('rec_output/rec_k_', num2str(max_reconstruction_degree),...
            '_', surface_name), f_rec, v_rec,'mode','ascii')
    catch
        stlWrite(strjoin('rec_output/rec_k_', num2str(max_reconstruction_degree),...
            '_', surface_name, ""), f_rec, v_rec,'mode','ascii')
    end
    if recursive_reconstruction
        %         This to export all the reconstruction degrees in different files
        for kk = 1:recursive_increments:max_reconstruction_degree-1
            v_rec_recrsive = real(rec_basis_func(:, 1:(kk+1)^2) * qm_k(1:(kk+1)^2, :));
            try
                stlWrite(strcat(['rec_output/rec_k_', num2str(kk),...
                    '_', surface_name]), f_rec, v_rec_recrsive,'mode','ascii')
            catch
                stlWrite(strjoin(['rec_output/rec_k_', num2str(kk),...
                    '_', surface_name], ""), f_rec, v_rec_recrsive,'mode','ascii')
            end
            fprintf("\nFile: (%s) has been saved.", ['rec_k_', num2str(kk),...
                '_', surface_name]);
        end
    end
    % This to export the final reconstruction degrees in one files
    try
        stlWrite(strcat(['rec_output/rec_k_', num2str(max_reconstruction_degree),...
            '_', surface_name]), v_rec(:, 1), v_rec(:, 2), v_rec(:, 3),'mode','ascii')
    catch
        stlWrite(strjoin(['rec_output/rec_k_', num2str(max_reconstruction_degree),...
            '_', surface_name], ""), v_rec(:, 1), v_rec(:, 2), v_rec(:, 3),'mode','ascii')
    end
    rec_height_map = v_rec(:, 3);
    if plot_figures
        paraview_patch(v_rec, f_rec, rec_height_map)
        axis equal
        title('Output surface')
        view([50 20])
    end
end
save(strcat(['rec_output/rec_k_', num2str(max_degree), '_', surface_name_wo_format,...
    '.mat']))
fprintf("\n\nFinished in %f min.\n", toc/60)
fprintf("\n\n Verts number: %3.0f", length(v(:, 1)))
fprintf("\n\n Equations number: %3.0f", length(v(:, 1)) * (max_degree+1)^2)
fprintf("\n\n (K+1)^2: %3.0f\n", (max_degree+1)^2)
fprintf("\n\n Total no. of coefficients 3(K+1)^2: %3.0f\n", 3*(max_degree+1)^2)
fprintf('The analysis ended at: %s\n', datestr(now,'HH:MM:SS.FFF'))