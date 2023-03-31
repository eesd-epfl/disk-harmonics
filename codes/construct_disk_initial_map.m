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

function initial_map = construct_disk_initial_map(v, f, flatten_method, aspect_ratio)

% Compute disk initial map using Tutte/locally authalic method

% for surface flattening, 1 = Tutte, 2 = locally authalic
if ~exist('flatten_method','var')
    flatten_method = 1;
end

% for the aspect ratio of the target shape
if ~exist('aspect_ratio','var')
    aspect_ratio = 1;
end


TR = TriRep(f,v); bdy_e = freeBoundary(TR); bdyid = bdy_e(:,1); % find bdy
bdy_length = sqrt(sum((v(bdyid,:) -  v(bdyid([2:end,1]),:)).^2,2));
partial_edge_sum = zeros(length(bdy_length),1);
% arc-length parameterization circular boundary constraint
for i = 2:length(bdy_length)
    for j = 1 : i-1
    partial_edge_sum(i) = partial_edge_sum(i) + bdy_length(j);
    end
end
theta = 2*pi.*partial_edge_sum/sum(bdy_length);
bdy = exp(theta*1i)';

bdy = real(bdy)*aspect_ratio + 1i*imag(bdy);

if flatten_method == 1
    % tutte map
  	initial_map = tutte_map(v,f,bdyid,bdy);
else
    % locally authalic map
    nv = length(v);
    M = authalic_matrix(v,f); 
    [mrow,mcol,mval] = find(M(bdyid,:));
    M = M - sparse(bdyid(mrow),mcol,mval,nv, nv) + ...
            sparse(bdyid,bdyid,ones(length(bdyid),1),nv, nv);
    c = zeros(nv,1); c(bdyid) = bdy;
    z = M \ c;
    initial_map = [real(z), imag(z)];
    
    % check bijectivity
    a = initial_map(f(:,2),:) - initial_map(f(:,1),:);
    b = initial_map(f(:,3),:) - initial_map(f(:,1),:);
    normal_direction = ((a(:,1).*b(:,2) - a(:,2).*b(:,1)) > 0);
    if length(unique(normal_direction)) > 1
       % contain overlaps, switch to tutte map
  	initial_map = tutte_map(v,f,bdyid,bdy);
    end
end
