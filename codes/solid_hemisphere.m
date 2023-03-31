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

function [v, f] = solid_hemisphere(N, a, b)
% Generate an STL mesh for a hemisphere. The hemisphere can be parbolic in
% shape, where a is the stretch along x and b is the stretch along y. N: is
% the number of refinement cycles.
% A scaled hemisphere can be produced as: [v, f] = SOLID_HEMISPHERE(N, a);
% A unit hemisphere can be produced as: [v, f] = SOLID_HEMISPHERE(N);
% A prabolic shape [v, f] = SOLID_HEMISPHERE(N, a, b);
if nargin == 0
    clc; clear; close all;
    N = 5;
    a = 1;
    b = 1;
elseif nargin == 1
    a = 1; b = 1;
elseif nargin == 2
    b = a;
end

% Generate a unit sphere N
[v, ~] = icosphere(N);
v(v(:, 3)>0, :) = [];

% Generate the top of the hemisphere
[v_d, ~] = uniform_disk_grid(0.1, N^2, 2); v_d(:, 3) = 0;

% Merge verticies
sz1 = size(v);
sz2 = size(v_d);
v(sz1(1)+1:sz1(1)+sz2(1), 1:3) = v_d;

% Scale the hemisphere if needed
v(:, 1) = a .* v(:, 1);
v(:, 2) = b .* v(:, 2);

% Remeshing to merge the surfaces
f = boundary(v);
[v,f] = clean_mesh(v,f,1);

if nargin == 0
    stlWrite("../tmp_hemisph.stl", f, v,'mode','ascii')
    scatter3(v(:,1), v(:,2), v(:,3))
    paraview_patch(v, f)
end
end