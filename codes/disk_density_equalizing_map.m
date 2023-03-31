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

function map = disk_density_equalizing_map(v,f)
% Compute an area-preserving parameterization of a simply-connected open
% surface onto a unit disk.
%
% v: nv x 3 vertex coordindates of the input mesh
% f: nf x 3 triangulations of the input mesh
% map: nv x 2 disk coordindates
%
% By Gary P. T. Choi, 2021

%% Density-equalizing map (Choi and Rycroft, SIIMS 2018)
bd = meshboundaries(f);
% initial map
initial_map = construct_disk_initial_map(v,f);
population = face_area(f,v);

% disk density-equalizing map
map = DEM_disk(v*sqrt(pi/sum(population)),f,population,initial_map);
% map = DEM_disk(v,f,population,initial_map);
map = map/max(abs(map(bd{1}))); % rescale to the unit disk

%% Check and enforce bijectivity 
mu = beltrami_coefficient2(initial_map,f,map);
mu(abs(mu)>0.8) = 0.8*mu(abs(mu)>0.8)./abs(mu(abs(mu)>0.8));
map = linear_beltrami_solver(initial_map,f,mu,bd{1},map(bd{1},:));

