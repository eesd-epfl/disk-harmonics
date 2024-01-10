function [H, K_G, I, II, norm_crv] = compute_curvatures(dD_rho, ddD_rho, dD_phi, ddD_phi, dD_rho_phi, Normals_n)
% Computes the coefficents of the first and second fundemental forms I and II. 
% I = [E|F|G]
% II = [e|f|g]
% Computes mean curvature H and Gaussian curvature K.
% Written by: Mahmoud Shaqfa @ 2024


% Compute I (first fundemental form)
I_E = sum(dD_rho.^2, 2);
I_F = sum(dD_rho .* dD_phi, 2);
I_G = sum(dD_phi.^2, 2);
I = [I_E, I_F, I_G];

% Compute II (second fundemental form)
II_e = sum(Normals_n .* ddD_rho, 2);
II_f = sum(Normals_n .* dD_rho_phi, 2);
II_g = sum(Normals_n .* ddD_phi, 2);
II = [II_e, II_f, II_g];

% Compute Gaussian curvature: det(dN)
K_G = (II_e .* II_g - II_f.^2);
norm_crv = (I_E .* I_G - I_F.^2);
norm_zero_idx = (norm_crv == 0);
norm_crv(norm_zero_idx) = 1;
K_G = K_G ./ norm_crv;
K_G(norm_zero_idx) = 0;

% Compute mean curvature: 1/2 * tr(dN)
H = 0.5 .* (II_e .* I_G - 2 .* II_f .* I_F + II_g .* I_E ); 
H = H ./ norm_crv;
H(norm_zero_idx) = 0;

end