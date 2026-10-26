function [U, eta, beta, tau0, tau1] = single_particle(ratio)
% SINGLE_PARTICLE  Exact single-particle constants for an elliptical surfer with
% aspect ratio a/b = ratio at unit area (ab = 1), surfactant profile g(s) = -cos(s),
% and s0 = pi/2: propulsion speed U, dipole strength eta, Bretherton constant beta, and
% the moments tau0, tau1 of the response to a uniform shift in the far-field
% concentration. See the exact solution and the far-field interactions in the manuscript.

s0 = pi/2;
a  = sqrt(ratio);  b = 1 ./ a;
lo = min(a, b);  m = 1 - (lo ./ max(a, b)).^2;
[K, E] = ellipke(m);

% k_i - k_0 along the propulsion axis x; a < b uses k_2 with the axes swapped
dk = -pi/8 * lo;
p = a > b;  dk(p) = lo(p) ./ (2*m(p)) .* (K(p) - E(p)) - lo(p)/2 .* K(p);
q = a < b;  dk(q) = lo(q) ./ (2*m(q)) .* (E(q) - (1 - m(q)) .* K(q)) - lo(q)/2 .* K(q);

U    = -(a + b) .* dk * (1 + cos(s0))^2 ./ (8 * a .* b);
eta  = pi * (a + b) * (1 + cos(s0))^2 / 4;
beta = (a.^2 - b.^2) ./ (a.^2 + b.^2);
tau0 = 2*pi ./ log(4 ./ (a + b));
tau1 = -(a + b) / 2 .* tau0;
end
