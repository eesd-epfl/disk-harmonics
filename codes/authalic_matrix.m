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

function L = authalic_matrix(v,f)

nv = length(v);
f1 = f(:,1); f2 = f(:,2); f3 = f(:,3);

% edge length
l = [sqrt(sum((v(f2,:) - v(f3,:)).^2,2)), sqrt(sum((v(f3,:) - v(f1,:)).^2,2)), sqrt(sum((v(f1,:) - v(f2,:)).^2,2))];
l1 = l(:,1); l2 = l(:,2); l3 = l(:,3);

% Heron's formula
s = (l1 + l2 + l3)*0.5;
area = sqrt(s.*(s-l1).*(s-l2).*(s-l3));

% weight
w12 = (l1.^2 + l2.^2 - l3.^2)./area/2;
w23 = (l2.^2 + l3.^2 - l1.^2)./area/2; 
w31 = (l1.^2 + l3.^2 - l2.^2)./area/2; 
diag1 = -w12./l2.^2 - w31./l3.^2;
diag2 = -w12./l1.^2 - w23./l3.^2;
diag3 = -w23./l2.^2 - w31./l1.^2;

% construct matrix
II = [f1; f2; f2; f3; f3; f1; f1; f2; f3];
JJ = [f3; f3; f1; f1; f2; f2; f1; f2; f3];
V = [w12./(l2.^2); w12./(l1.^2); w23./(l3.^2); w23./(l2.^2); w31./(l1.^2); w31./(l3.^2); diag1; diag2; diag3];
L = sparse(II,JJ,V,nv,nv);

end