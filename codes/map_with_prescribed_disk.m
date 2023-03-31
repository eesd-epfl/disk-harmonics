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
% Author of this file: Gary CHoi

function map = map_with_prescribed_disk(v,f,vid_patch)

% Compute a planar parameterization of a genus-0 closed surface, with
% a prescribed open surface patch of it mapped to the unit disk
%
% Input:
% v: nv x 3 vertex coordinates of a genus-0 closed triangle mesh
% f: nf x 3 triangulations of a genus-0 closed triangle mesh
% vid_patch: the vertex indices of the open surface patch to be considered
% 
% Output:
% map: nv x 2 vertex coordinates of the planar parameterization
%
% Written by Gary Choi, 2021

%% For test purposes
if nargin == 0
    %%
    close all; clc;
   % Assuming that all the sub-directories are loaded
   [v,f] = stlRead("cube_refined_sample_2.stl");
   filter_query = strcat('vid_patch = find(v(:,3)>0.577042-0.01);'); % for the donor mesh (scaled cube)
   eval(filter_query);
end

%% extract the open patch
% find faces for which all three vertices are in the open patch
fid_patch = find(ismember(f(:,1),vid_patch).*ismember(f(:,2),vid_patch).*ismember(f(:,3),vid_patch));
[f_patch1,v_patch1,father1] = gpp_clean_mesh(f(fid_patch,:),v);
[f_patch2,v_patch2,father2] = gpp_clean_mesh(f(setdiff(1:length(f),fid_patch),:),v);

bdyid_patch = meshboundaries(f_patch1); bdyid_patch = bdyid_patch{1};
bdyid_patch2 = meshboundaries(f_patch2); bdyid_patch2 = bdyid_patch2{1};

% build the boundary correspondence
bdy_correspondence = zeros(length(bdyid_patch2),1);
for i = 1:length(bdyid_patch2)
    [~,id] = min(abs(v_patch1(bdyid_patch,1) - v_patch2(bdyid_patch2(i),1)) + ...
                 abs(v_patch1(bdyid_patch,2) - v_patch2(bdyid_patch2(i),2)) + ...
                 abs(v_patch1(bdyid_patch,3) - v_patch2(bdyid_patch2(i),3)));
    bdy_correspondence(i) = id;
end

if nargin == 0
    plot_mesh(v_patch1,f_patch1); 
end
%% map the open patch onto the unit disk
map_patch1 = disk_density_equalizing_map(v_patch1,f_patch1);
fz1 = complex(map_patch1(:,1),map_patch1(:,2));
if nargin == 0
    plot_mesh(map_patch1,f_patch1);
end
%% map the other patch to the exterior of the disk
disk_bdy = fz1(bdyid_patch);

% disk harmonic map
M = cotangent_laplacian(v_patch2,f_patch2);
[mrow,mcol,mval] = find(M(bdyid_patch2,:));
M = M - sparse(bdyid_patch2(mrow),mcol,mval,length(v_patch2), length(v_patch2)) + ...
        sparse(bdyid_patch2,bdyid_patch2,ones(length(bdyid_patch2),1),length(v_patch2), length(v_patch2));
c = zeros(length(v_patch2),1); 
c(bdyid_patch2) = disk_bdy(bdy_correspondence);
fz2 = 1./conj(M \ c); % apply a reflection to map everything to the exterior of the disk
if nargin == 0
    plot_mesh([real(fz2),imag(fz2)],f_patch2); 
end
%% combine the results and obtain the final planar map
fz = zeros(length(v),1);
fz(father1) = fz1;
fz(father2) = fz2;
map = [real(fz),imag(fz)];
if nargin == 0
    plot_mesh(map,f);
    axis on
    xlabel('x');
    ylabel('y');
end
end

