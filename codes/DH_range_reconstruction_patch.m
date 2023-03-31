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

function recursive_reconstruction = DH_range_reconstruction_patch(k_range, qm_k, eigen_table, rhos, phis)
% Icosahedron dome reconstruction
k_min = min(k_range);
k_max = max(k_range);

C_mat = zeros(length(rhos), (k_max+1)^2);
reconstruction = zeros(length(rhos), 3);
recursive_reconstruction = zeros(length(rhos), 3, k_max + 1);
for n = k_min:k_max
    for l = 0:n
        for m = -l:1:l
            lmn = eigen_table(l+1, abs(m)+1);
            Normalization = 0.5 * (1 - (abs(m)/lmn)^2) * besselj(abs(m), lmn)^2;
            Normalization = (1/sqrt(Normalization)) * (1/sqrt(2*pi));
            C_mat(:, l^2 + l + m + 1) = Normalization .* besselj(abs(m), lmn*rhos) .* exp(1i*abs(m).*phis);
            if m < 0
                C_mat(:, l^2 + l + m + 1) = conj(C_mat(:, l^2 + l + m + 1)) .* (-1)^m;
            end
        end
    end
    
    for ii = n^2 + 1:n^2 + 2*n + 1
        reconstruction(:, 1) = reconstruction(:, 1) + C_mat(:, ii) .* qm_k(ii, 1);
        reconstruction(:, 2) = reconstruction(:, 2) + C_mat(:, ii) .* qm_k(ii, 2);
        reconstruction(:, 3) = reconstruction(:, 3) + C_mat(:, ii) .* qm_k(ii, 3);
    end
    recursive_reconstruction(:, :, n+1) = reconstruction;
    fprintf("\nCompleted: %3.2f %% for n = %3.0f", ((n+1)^2/(k_max + 1)^2*100), n)
end
end