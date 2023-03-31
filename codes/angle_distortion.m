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

function distortion = angle_distortion(v,f,map)

% Calculate and visualize the angle difference (angle_map - angle_v)
% 
% Input:
% v: nv x (1/2/3) vertex coordinates 
% f: nf x 3 triangulations 
% map: nv x (1/2/3) vertex coordinates of the mapping result
%
% Output:
% distortion: 3*nf x 1 angle differences
% 
% If you use this code in your own work, please cite the following paper:
% [1] P. T. Choi and L. M. Lui, 
%     "Fast Disk Conformal Parameterization of Simply-Connected Open Surfaces."
%     Journal of Scientific Computing, 65(3), pp. 1065-1090, 2015.
%
% Copyright (c) 2014-2018, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

nv = length(v);
nv2 = length(map);

if nv ~= nv2
    error('Error: The two meshes are of different size.');
end

if size(v,2) == 1
    v = [real(v),imag(v),zeros(length(v),1)];
elseif size(v,2) == 2
    v = [v,zeros(length(v),1)];
end

if size(map,2) == 1
    map = [real(map),imag(map),zeros(length(map),1)];
elseif size(map,2) == 2
    map = [map,zeros(length(map),1)];
end

f1 = f(:,1); f2 = f(:,2); f3 = f(:,3);

% calculate angles on v
a3=[v(f1,1)-v(f3,1), v(f1,2)-v(f3,2), v(f1,3)-v(f3,3)];
b3=[v(f2,1)-v(f3,1), v(f2,2)-v(f3,2), v(f2,3)-v(f3,3)];
a1=[v(f2,1)-v(f1,1), v(f2,2)-v(f1,2), v(f2,3)-v(f1,3)];
b1=[v(f3,1)-v(f1,1), v(f3,2)-v(f1,2), v(f3,3)-v(f1,3)];
a2=[v(f3,1)-v(f2,1), v(f3,2)-v(f2,2), v(f3,3)-v(f2,3)];
b2=[v(f1,1)-v(f2,1), v(f1,2)-v(f2,2), v(f1,3)-v(f2,3)];
vcos1=(a1(:,1).*b1(:,1)+a1(:,2).*b1(:,2)+a1(:,3).*b1(:,3))./ ...
    (sqrt(a1(:,1).^2+a1(:,2).^2+a1(:,3).^2).*sqrt(b1(:,1).^2+b1(:,2).^2+b1(:,3).^2));
vcos2=(a2(:,1).*b2(:,1)+a2(:,2).*b2(:,2)+a2(:,3).*b2(:,3))./...
    (sqrt(a2(:,1).^2+a2(:,2).^2+a2(:,3).^2).*sqrt(b2(:,1).^2+b2(:,2).^2+b2(:,3).^2));
vcos3=(a3(:,1).*b3(:,1)+a3(:,2).*b3(:,2)+a3(:,3).*b3(:,3))./...
    (sqrt(a3(:,1).^2+a3(:,2).^2+a3(:,3).^2).*sqrt(b3(:,1).^2+b3(:,2).^2+b3(:,3).^2));
    
% calculate angles on map
c3=[map(f1,1)-map(f3,1), map(f1,2)-map(f3,2), map(f1,3)-map(f3,3)];
d3=[map(f2,1)-map(f3,1), map(f2,2)-map(f3,2), map(f2,3)-map(f3,3)];
c1=[map(f2,1)-map(f1,1), map(f2,2)-map(f1,2), map(f2,3)-map(f1,3)];
d1=[map(f3,1)-map(f1,1), map(f3,2)-map(f1,2), map(f3,3)-map(f1,3)];
c2=[map(f3,1)-map(f2,1), map(f3,2)-map(f2,2), map(f3,3)-map(f2,3)];
d2=[map(f1,1)-map(f2,1), map(f1,2)-map(f2,2), map(f1,3)-map(f2,3)];
mapcos1=(c1(:,1).*d1(:,1)+c1(:,2).*d1(:,2)+c1(:,3).*d1(:,3))./...
    (sqrt(c1(:,1).^2+c1(:,2).^2+c1(:,3).^2).*sqrt(d1(:,1).^2+d1(:,2).^2+d1(:,3).^2));
mapcos2=(c2(:,1).*d2(:,1)+c2(:,2).*d2(:,2)+c2(:,3).*d2(:,3))./...
    (sqrt(c2(:,1).^2+c2(:,2).^2+c2(:,3).^2).*sqrt(d2(:,1).^2+d2(:,2).^2+d2(:,3).^2));
mapcos3=(c3(:,1).*d3(:,1)+c3(:,2).*d3(:,2)+c3(:,3).*d3(:,3))./...
    (sqrt(c3(:,1).^2+c3(:,2).^2+c3(:,3).^2).*sqrt(d3(:,1).^2+d3(:,2).^2+d3(:,3).^2));
    
% calculate the angle difference
distortion = rad2deg(mean(abs(acos([mapcos1;mapcos2;mapcos3])-acos([vcos1;vcos2;vcos3]))));
