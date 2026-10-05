% vonmises_dispersion.m — solve the dispersion relation for the von Mises state
%   1 = (i eta/2) \int (khat.p) h(theta) psi_0 / (sigma + i U k (khat.p)) dtheta,
%   h = -(1+beta) + 2 beta (khat.p)^2 + zeta beta sin(theta) (khat.p_perp)(khat.p),
% for the unstable growth rates sigma with Newton's method, and plot the maximal
% growth rate in the k-phi plane. Initial guesses are the unstable eigenvalues
% of the discretized operator (as in vonmises_stability.m).

close all

fig = journal_figure([6.75 2.5], 2, 'tex');
a = 2; b = 1/a;
beta = (a^2 - b^2)/(a^2 + b^2);
U = 1;
eta = 1;

N = 64;
zeta = 15;
kh = linspace(0, 5, 32);
ph = linspace(0, pi/2, 32);
[K, PH] = meshgrid(kh, ph);
Smax = arrayfun(@(k, p) max_growth(k, p, zeta, eta, U, beta, N), K, PH);

ax = gca; hold(ax, 'on');
contourf(ax, K, PH, Smax, 8); axis(ax, 'xy', 'tight'); colorbar(ax);
pc = linspace(1e-3, pi/2, 200);
plot(ax, eta*(1 - beta*sin(pc).^2).^2 ./ (8*beta*U*sin(pc).^2), pc, 'k--', 'LineWidth', 1.5, 'DisplayName', '$k_c$ ($\zeta \rightarrow \infty$)');

xlim(ax, kh([1 end]));
ylim(ax, ph([1 end]));

yticks(ax, [0 pi/8 pi/4 3*pi/8 pi/2]);
yticklabels(ax, {'0', '\pi/8', '\pi/4', '3\pi/8', '\pi/2'});

xlabel(ax, '\it{k}');
ylabel(ax, '\phi');

cmap = cmocean('matter');
cmap = [1 1 1; cmap];
colormap(ax, cmap);

function smax = max_growth(k, phi, zeta, eta, U, beta, N)

    sigma_0 = operator_eigs(k, phi, zeta, eta, U, beta, N);
    sigma_0 = sigma_0(real(sigma_0) > 1e-8);
    f = @(s) dispersion(s, k, phi, zeta, eta, U, beta);
    sigma = arrayfun(@(s) newton(f, s), sigma_0);
    residual = arrayfun(@(s) abs(f(s)), sigma);
    smax = max([0; real(sigma(residual < 1e-8 & real(sigma) > 0))]);

end

function [D, dD] = dispersion(sigma, k, phi, zeta, eta, U, beta)

    c = @(t) cos(t - phi);      % khat . p
    s = @(t) -sin(t - phi);     % khat . p_perp
    psi_0 = @(t) exp(zeta*(cos(t) - 1)) / (2*pi*besseli(0, zeta, 1));
    h = @(t) -(1 + beta) + 2*beta*c(t).^2 + zeta*beta*sin(t).*s(t).*c(t);
    den = @(t) sigma + 1i*U*k*c(t);

    D = 1 - (1i*eta/2) * integral(@(t) c(t).*h(t).*psi_0(t)./den(t), 0, 2*pi, 'RelTol', 1e-8);
    dD = (1i*eta/2) * integral(@(t) c(t).*h(t).*psi_0(t)./den(t).^2, 0, 2*pi, 'RelTol', 1e-8);

end

function s = newton(f, s)

    for it = 1 : 50
        [D, dD] = f(s);
        ds = D / dD;
        s = s - ds;
        if abs(ds) < 1e-8 * max(1, abs(s)), break, end
    end

end

function sigma = operator_eigs(k, phi, zeta, eta, U, beta, N)

    t = linspace(-pi, pi, N+1)';
    t = t(1 : end-1);
    w = t(2) - t(1);
    p = [cos(t), sin(t)];
    p_perp = [-sin(t), cos(t)];
    khat = [cos(phi); sin(phi)];

    psi_0 = exp(zeta*(cos(t) - 1)) / (2*pi*besseli(0, zeta, 1));
    dpsi_0 = -zeta*sin(t) .* psi_0;

    m_1 = w * p.';
    ku_1 = (eta/2) * khat * (khat.' * m_1);
    E_1 = @(a, b) (1i/2) * ((a*khat) .* (b*ku_1) + (a*ku_1) .* (b*khat));

    tdot_1 = beta * E_1(p_perp, p);
    dtdot_1 = beta * (E_1(p_perp, p_perp) - E_1(p, p));

    L = -1i*k*diag(U*p*khat) - psi_0 * (1i*khat.' * ku_1) - (dtdot_1 .* psi_0 + tdot_1 .* dpsi_0);
    sigma = eig(L);

end
