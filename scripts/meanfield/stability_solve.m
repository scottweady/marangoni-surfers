% stability_solve.m — maximal growth rate of the linearized mean-field model
% about the von Mises state psi_0 ∝ exp(kappa cos theta), as a function of
% wavenumber k and wavevector angle phi, for kappa = 5, 10, 15. Growth rates
% are roots of the dispersion relation (see manuscript/vonmises_dispersion.m).
% Saves data/stability/growth_rate.mat, read by figs/fig_stability.m.

warning('off')
a = 2; b = 1/a;
lambda = (a^2 - b^2)/(a^2 + b^2);
U = 1;
beta = 1;

N = 32;
kappas = [5 10 15];
kh = linspace(0, 6, 200);
ph = linspace(0, pi/2, 200);
Smax = cell(1, numel(kappas));
tic
for j = 1:numel(kappas)
    Smax{j} = zeros(numel(ph), numel(kh));
    for i = 1:numel(ph)
        Smax{j}(i,:) = growth_row(kh, ph(i), kappas(j), beta, U, lambda, N);
        nk = find(Smax{j}(i,:) > 0, 1, 'last');
        fprintf('kappa = %g  phi %3d/%d (%.3f)  unstable to k = %.2f  max sigma = %.3g  [%.0f s]\n', ...
            kappas(j), i, numel(ph), ph(i), max([0 kh(nk)]), max(Smax{j}(i,:)), toc);
    end
end

outDir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data', 'stability');
if ~exist(outDir, 'dir'), mkdir(outDir); end
save(fullfile(outDir, 'growth_rate.mat'), 'kh', 'ph', 'kappas', 'Smax', 'a', 'lambda', 'U', 'beta', 'N');
fprintf('Saved %s\n', fullfile(outDir, 'growth_rate.mat'));

contourf(kh,ph,Smax{1})
function srow = growth_row(kh, phi, kappa, beta, U, lambda, N)
% Max growth rate along a row of fixed phi. Marches up in k, continuing the
% root from the previous k; falls back to the full search when that fails, and
% stops once nstop consecutive points are stable after an unstable band.

    tol   = 1e-6;   % growth rate below which a point doesn't mark the row unstable
    nstop = 10;
    srow = zeros(size(kh));
    s = NaN;
    unstable = false;
    nstable  = 0;
    for n = 1:numel(kh)
        ok = false;
        if ~isnan(s)
            f = @(s) dispersion(s, kh(n), phi, kappa, beta, U, lambda);
            [s, ok] = newton(f, s);
        end
        if ~ok, s = max_growth(kh(n), phi, kappa, beta, U, lambda, N); end
        if isnan(s)
            nstable = nstable + 1;
            if unstable && nstable >= nstop, break, end
            continue
        end
        nstable = 0;
        unstable = unstable || real(s) > tol;
        srow(n) = real(s);
    end

end

function smax = max_growth(k, phi, kappa, beta, U, lambda, N)
% Root with the largest real part (NaN if none with Re > 0), refined by Newton
% from the eigenvalues of the discretized operator.

    sigma_0 = operator_eigs(k, phi, kappa, beta, U, lambda, N);
    sigma_0 = sigma_0(real(sigma_0) > 1e-12);
    f = @(s) dispersion(s, k, phi, kappa, beta, U, lambda);
    smax = NaN;
    for s0 = sigma_0.'
        [s, ok] = newton(f, s0);
        if ok && (isnan(smax) || real(s) > real(smax)), smax = s; end
    end

end

function [D, dD] = dispersion(sigma, k, phi, kappa, beta, U, lambda)

    c = @(t) cos(t - phi);      % khat . p
    s = @(t) -sin(t - phi);     % khat . p_perp
    psi_0 = @(t) exp(kappa*(cos(t) - 1)) / (2*pi*besseli(0, kappa, 1));
    h = @(t) -(1 + lambda) + 2*lambda*c(t).^2 + kappa*lambda*sin(t).*s(t).*c(t);
    den = @(t) sigma + 1i*U*k*c(t);

    D = 1 - (1i*beta/2) * integral(@(t) c(t).*h(t).*psi_0(t)./den(t), 0, 2*pi, 'RelTol', 1e-6);
    dD = (1i*beta/2) * integral(@(t) c(t).*h(t).*psi_0(t)./den(t).^2, 0, 2*pi, 'RelTol', 1e-6);

end

function [s, ok] = newton(f, s)
% Gives up once the iterate leaves Re s > 0, since only those roots are kept.

    ok = false;
    for it = 1 : 25
        [D, dD] = f(s);
        ds = D / dD;
        s = s - ds;
        if real(s) <= 0, return, end
        if abs(ds) < 1e-12 * max(1, abs(s)), ok = true; return, end
    end

    fprintf('Newton did not converge for initial guess %g\n', s);

end

function sigma = operator_eigs(k, phi, kappa, beta, U, lambda, N)

    t = linspace(-pi, pi, N+1)';
    t = t(1 : end-1);
    w = t(2) - t(1);
    p = [cos(t), sin(t)];
    p_perp = [-sin(t), cos(t)];
    khat = [cos(phi); sin(phi)];

    psi_0 = exp(kappa*(cos(t) - 1)) / (2*pi*besseli(0, kappa, 1));
    dpsi_0 = -kappa*sin(t) .* psi_0;

    m_1 = w * p.';
    ku_1 = (beta/2) * khat * (khat.' * m_1);
    E_1 = @(a, b) (1i/2) * ((a*khat) .* (b*ku_1) + (a*ku_1) .* (b*khat));

    tdot_1 = lambda * E_1(p_perp, p);
    dtdot_1 = lambda * (E_1(p_perp, p_perp) - E_1(p, p));

    L = -1i*k*diag(U*p*khat) - psi_0 * (1i*khat.' * ku_1) - (dtdot_1 .* psi_0 + tdot_1 .* dpsi_0);
    sigma = eig(L);

end
