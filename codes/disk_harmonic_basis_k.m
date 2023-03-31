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

function SCH_basis = disk_harmonic_basis_k(K, eigen_table, rhos, phis)
% Imeplementation: Mahmoud Shaqfa, Gary Choi, Katrin Beyer.
% It returns a table for positive and negative orders for a certain index K
% The number of the columns in the matrix is 2K+1 for (-m:1:m) orders
% Note: rhos and phis matrix msut be identical in size
% All the angles here must be in radians.
% For our applications the following variables shall be fixed:
% rhos and phis must be square matrix or column matrix (row matrix will give you an error!)
BC = "even";
normalization_method = "Hainse";

% Initializing and filling up
ms = 0:K;

% Check inputs' size (1D and 2D matrix)
sz = length(size(rhos)); % For reconstruction matricies only
if isempty(find(size(rhos) == 1, 1))
    sz = sz + 1;
end
if sz == 2 % For normal analysis or reconstruction using geodesic domes
    SCH_basis = zeros(length(rhos), (2*K+1));
    positive_ms = repmat(ms, length(rhos), 1);
    % Compute and fill up the postive orders
    SCH_basis(:, K+1:end) = fractional_associated_Legendre_functions(theta_c,...
        K, BC, eigen_table, rhos, normalization_method, N_eps) .* exp(1i .* positive_ms .* repmat(phis, 1, K+1));
    % Compute and fill up the negative orders
    SCH_basis(:, 1:K) = flip(conj(SCH_basis(:, K+2:end)) .* (-1).^ positive_ms(:, 2:end), 2);
% elseif sz == 3 % For the meshgrid reconstruction case only or 3D render of the basis
%     temp_sz = size(rhos); temp_sz(3) = (2*K+1);
%     SCH_basis = zeros(temp_sz);
%     % Prepare the ms
%     positive_ms = zeros(temp_sz(1), temp_sz(2), K+1);
%     for ii = ms
%         positive_ms(:, :, ii+1) = ones(temp_sz(1), temp_sz(2)) .* ii;
%     end
%     % Compute and fill up the postive orders
%     SCH_basis(:, :, K+1:end) = fractional_associated_Legendre_functions(deg2rad(theta_c),...
%         K, BC, eigen_table, rhos, normalization_method, N_eps) .* exp(1i .* positive_ms .* phis);
%     % Compute and fill up the negative orders
%     SCH_basis(:, :, 1:K) = flip(conj(SCH_basis(:, :, K+2:end)) .* (-1).^ positive_ms(:, :, 2:end), 3);
end
end