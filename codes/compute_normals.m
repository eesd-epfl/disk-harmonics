function [Normals, Normals_n] = compute_normals(dD_phi, dD_rho)
% Compute analytic normals
% Mahmoud Shaqfa @ 2024

Normals = cross(dD_phi, dD_rho, 2);
norm_ = sqrt(sum(Normals.^2, 2));
Normals_n = Normals ./ norm_;
end