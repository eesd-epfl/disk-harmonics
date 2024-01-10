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
% Author of this file: Mahmoud Shaqfa @ 2024

function [dD_rho, ddD_rho, dD_phi, ddD_phi, dD_rho_phi] = derivatives_disk_harmonic(K, qm_k, eigen_table, rhos, phis)
% Imeplementation: Mahmoud Shaqfa
% Derivatives of Bessel-Fourier basis functions for the disk harmonics (\rho, \phi).
% Note: rhos and phis matrix msut be identical in size
% All the angles here must be in radians.
% For our applications the following variables shall be fixed:
% rhos and phis must be square matrix or column matrix (row matrix will give you an error!)
% BC = "Neumann";

dD_rho = zeros(size(rhos, 1), 3);
ddD_rho = zeros(size(rhos, 1), 3);
dD_phi = zeros(size(rhos, 1), 3);
ddD_phi = zeros(size(rhos, 1), 3);
dD_rho_phi = zeros(size(rhos, 1), 3);

for n = 0:K
    %   For positive orders
    for m = -n:n
        lmn = eigen_table(n+1, abs(m)+1);
        idx = n^2 + n + m + 1;
        % Neumann BCs
        Normalization = 0.5 * (1 - (abs(m)/lmn)^2) * besselj(abs(m), lmn)^2;
        Normalization = 1/sqrt(Normalization);          % Bessel normalization (for Neumann BC only)
        Normalization = Normalization * (1/sqrt(2*pi)); % Add Fourier normalization
        [dy_Jn, ddy_Jn] = derivative_bessel_basis(abs(m), lmn, rhos);

        if m >= 0
            % Partial derivatives with respect to \rho
            dD_rho(:, 1) = dD_rho(:, 1) + Normalization .* dy_Jn .* exp(1i*abs(m).*phis) .* qm_k(idx, 1); %x
            dD_rho(:, 2) = dD_rho(:, 2) + Normalization .* dy_Jn .* exp(1i*abs(m).*phis) .* qm_k(idx, 2); %y
            dD_rho(:, 3) = dD_rho(:, 3) + Normalization .* dy_Jn .* exp(1i*abs(m).*phis) .* qm_k(idx, 3); %z

            ddD_rho(:, 1) = ddD_rho(:, 1) + Normalization .* ddy_Jn .* exp(1i*abs(m).*phis) .* qm_k(idx, 1); %x
            ddD_rho(:, 2) = ddD_rho(:, 2) + Normalization .* ddy_Jn .* exp(1i*abs(m).*phis) .* qm_k(idx, 2); %y
            ddD_rho(:, 3) = ddD_rho(:, 3) + Normalization .* ddy_Jn .* exp(1i*abs(m).*phis) .* qm_k(idx, 3); %z

            % Partial derivatives with respect to \phi
            dD_phi(:, 1) = dD_phi(:, 1) + 1i * abs(m) * Normalization .* besselj(abs(m), lmn .* rhos) .* exp(1i*abs(m).*phis) .* qm_k(idx, 1); %x
            dD_phi(:, 2) = dD_phi(:, 2) + 1i * abs(m) * Normalization .* besselj(abs(m), lmn .* rhos) .* exp(1i*abs(m).*phis) .* qm_k(idx, 2); %y
            dD_phi(:, 3) = dD_phi(:, 3) + 1i * abs(m) * Normalization .* besselj(abs(m), lmn .* rhos) .* exp(1i*abs(m).*phis) .* qm_k(idx, 3); %z

            ddD_phi(:, 1) = ddD_phi(:, 1) + (-1) * m^2 * Normalization .* besselj(abs(m), lmn .* rhos) .* exp(1i*abs(m).*phis) .* qm_k(idx, 1); %x
            ddD_phi(:, 2) = ddD_phi(:, 2) + (-1) * m^2 * Normalization .* besselj(abs(m), lmn .* rhos) .* exp(1i*abs(m).*phis) .* qm_k(idx, 2); %y
            ddD_phi(:, 3) = ddD_phi(:, 3) + (-1) * m^2 * Normalization .* besselj(abs(m), lmn .* rhos) .* exp(1i*abs(m).*phis) .* qm_k(idx, 3); %z

            % Mixed Partial derivatives with respect to \phi and \rho
            dD_rho_phi(:, 1) =  dD_rho_phi(:, 1) + 1i * abs(m) * Normalization .* dy_Jn .* exp(1i*abs(m).*phis) .* qm_k(idx, 1); %x
            dD_rho_phi(:, 2) =  dD_rho_phi(:, 2) + 1i * abs(m) * Normalization .* dy_Jn .* exp(1i*abs(m).*phis) .* qm_k(idx, 2); %y
            dD_rho_phi(:, 3) =  dD_rho_phi(:, 3) + 1i * abs(m) * Normalization .* dy_Jn .* exp(1i*abs(m).*phis) .* qm_k(idx, 3); %z
        elseif m < 0
            % Partial derivatives with respect to \rho
            dD_rho(:, 1) = dD_rho(:, 1) + (-1)^abs(m) * Normalization .* dy_Jn .* exp(-1i*abs(m).*phis) .* qm_k(idx, 1); %x
            dD_rho(:, 2) = dD_rho(:, 2) + (-1)^abs(m) * Normalization .* dy_Jn .* exp(-1i*abs(m).*phis) .* qm_k(idx, 2); %y
            dD_rho(:, 3) = dD_rho(:, 3) + (-1)^abs(m) * Normalization .* dy_Jn .* exp(-1i*abs(m).*phis) .* qm_k(idx, 3); %z
            
            ddD_rho(:, 1) = ddD_rho(:, 1) + (-1)^abs(m) * Normalization .* ddy_Jn .* exp(-1i*abs(m).*phis) .* qm_k(idx, 1); %x
            ddD_rho(:, 2) = ddD_rho(:, 2) + (-1)^abs(m) * Normalization .* ddy_Jn .* exp(-1i*abs(m).*phis) .* qm_k(idx, 2); %y
            ddD_rho(:, 3) = ddD_rho(:, 3) + (-1)^abs(m) * Normalization .* ddy_Jn .* exp(-1i*abs(m).*phis) .* qm_k(idx, 3); %z
            
            % Partial derivatives with respect to \phi
            dD_phi(:, 1) = dD_phi(:, 1) + (-1)^abs(m) * -1i * abs(m) * Normalization .* besselj(abs(m), lmn .* rhos) .* exp(-1i*abs(m).*phis) .* qm_k(idx, 1); %x
            dD_phi(:, 2) = dD_phi(:, 2) + (-1)^abs(m) * -1i * abs(m) * Normalization .* besselj(abs(m), lmn .* rhos) .* exp(-1i*abs(m).*phis) .* qm_k(idx, 2); %y
            dD_phi(:, 3) = dD_phi(:, 3) + (-1)^abs(m) * -1i * abs(m) * Normalization .* besselj(abs(m), lmn .* rhos) .* exp(-1i*abs(m).*phis) .* qm_k(idx, 3); %z
            
            ddD_phi(:, 1) = ddD_phi(:, 1) + (-1)^abs(m+1) * (-1) * m^2 * Normalization .* besselj(abs(m), lmn .* rhos) .* exp(-1i*abs(m).*phis) .* qm_k(idx, 1); %x
            ddD_phi(:, 2) = ddD_phi(:, 2) + (-1)^abs(m+1) * (-1) * m^2 * Normalization .* besselj(abs(m), lmn .* rhos) .* exp(-1i*abs(m).*phis) .* qm_k(idx, 2); %y
            ddD_phi(:, 3) = ddD_phi(:, 3) + (-1)^abs(m+1) * (-1) * m^2 * Normalization .* besselj(abs(m), lmn .* rhos) .* exp(-1i*abs(m).*phis) .* qm_k(idx, 3); %z
            
            % Mixed Partial derivatives with respect to \phi and \rho
            dD_rho_phi(:, 1) =  dD_rho_phi(:, 1) + (-1)^abs(m) * -1i * abs(m) * Normalization .* dy_Jn .* exp(-1i*abs(m).*phis) .* qm_k(idx, 1); %x
            dD_rho_phi(:, 2) =  dD_rho_phi(:, 2) + (-1)^abs(m) * -1i * abs(m) * Normalization .* dy_Jn .* exp(-1i*abs(m).*phis) .* qm_k(idx, 2); %y
            dD_rho_phi(:, 3) =  dD_rho_phi(:, 3) + (-1)^abs(m) * -1i * abs(m) * Normalization .* dy_Jn .* exp(-1i*abs(m).*phis) .* qm_k(idx, 3); %z
        end
    end
    fprintf("\nCompleted: %3.2f %% for n = %3.0f", ((n+1)^2/(K+1)^2*100), n)
end