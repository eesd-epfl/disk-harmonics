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
% Author of this file: Mahmoud Shaqfa

function eigen_table = Sturm_Liouville_eigenvalues_v2(K)
% Authors: Mahmoud Shaqfa, Gary Choi, Katrin Beyer.
% This is a solver to return the first n-roots of Bessel polynomials
% by applying either Dirichlet or Neumann Boundary conditions.
% K: (Postive integer) the index of the m orders.
% BC: the boundary conditions type. For Dirichlet we apply constant
% boundary of the expanded function (odd basis). For Neumann we apply
% boundaries for the first derivative and equate it by zero (even basis).
% BC is "neumann" or "dirichlet".
% The function returns KxM table of roots.
% example:
%         eigen_table = Sturm_Liouville_eigenvalues(15, 'neumann')
% solver_increment_diag: the diagonal solver for all m = k; this should be
% smaller than the solver_increment_tri for the rest of the m > 0 roots.
if nargin == 0
    close all; clc; tic
    K = 70;
end
eigen_table = zeros(length(0:1:K));
%% Start the solver here
% Find the first K-roots (eigenvalues) at \theta_o

% Solve for diagonals
% Find the first non-zero roots here.
for k = 1:K+1
    eigen_table(k:end, k) = besselzeroj(k-1, K+2-k, 'd');
end

if nargin == 0
    toc
    save(strcat("eigenvalues_k_", num2str(K), ".mat"), 'eigen_table')
    fprintf("\nThe solver time: %1.5f min\n", toc/60)
    fig = figure;
%     set(fig,'renderer','painters')
    b = bar3(eigen_table'); colorbar; xlabel("k"); ylabel("m"); zlabel("l(m)")
    for k = 1:length(b)
        zdata = b(k).ZData;
        b(k).CData = zdata;
        b(k).FaceColor = 'interp';
    end
    set(gca,'XTick', 0:10:K);
    set(gca,'YTick', 0:10:K);
    set(gca,'fontname','Amiri')  % Set it to times
    set(gca,'FontSize',15)
    view(2)
end
end