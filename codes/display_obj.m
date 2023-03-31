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

function display_obj(obj,texture)
%
% function display_obj(obj,texture)
%
% displays an obj structure with texture
%
% INPUTS:  obj:     object data
%                   - obj.v:    vertices
%                   - obj.vt:   texture coordinates
%                   - obj.f.v:  face definition vertices
%                   - obj.f.vt: face definition texture
%
%       : texture -  texture image full path
%
% Author: Bernard Abayowa
% University of Dayton
% 6/16/08
texture = imread(texture);
texture_img = flipdim(texture,1);
[sy sx sz] = size(texture_img);
texture_img =  reshape(texture_img,sy*sx,sz);
% make image 3D if grayscale
if sz == 1
    texture_img = repmat(texture_img,1,3);
end
% select what texture correspond to each vertex according to face
% definition
[vertex_idx fv_idx] = unique(obj.f.v);
texture_idx = obj.f.vt(fv_idx);
x = abs(round(obj.vt(:,1)*(sx-1)))+1;
y = abs(round(obj.vt(:,2)*(sy-1)))+1;
xy = sub2ind([sy sx],y,x);
texture_pts = xy(texture_idx);
tval = double(texture_img(texture_pts,:))/255;
% display object
figure; patch('vertices',obj.v,'faces',obj.f.v,'FaceVertexCData', tval);
shading interp
colormap gray(256);
lighting phong;
camproj('perspective');
axis square; 
axis off;
axis equal
axis tight;
cameramenu