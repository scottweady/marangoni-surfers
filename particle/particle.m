
close all

dt = 0.1;
tf = 1000.0;

a = 0.25; b = 1.0/a; beta = 0.048896; uS = 0.003081;
% a = 1.0;  b = 1.0;  beta = 1.570796;  uS = 0.098175;
% a = 2.0;  b = 0.5;  beta = 1.689178;  uS = 0.128646;

lam = (a^2 - b^2) / (a^2 + b^2);
R0  = beta / (4*pi * uS) * (1 - lam/2);

R_vals  = R0 * linspace(0.5, 1.5, 7);
R_vals = [R0];
th_vals = linspace(-pi, pi, 10);
% th_vals = [pi/2];
cmap    = lines(numel(th_vals));

figure; hold on; box on;
for k = 1:numel(R_vals)
    R = R_vals(k);
    x = [-R/2, 0.0; R/2, 0.0];
    for j = 1:numel(th_vals)
        th  = [-(pi/2+0.5*th_vals(j)); pi/2+0.5*th_vals(j)];
        vis = 'off';
        if k == 1; vis = 'on'; end
        [t_hist, sep_hist] = simulate(x, th, beta, lam, uS, dt, tf);
        plot(t_hist, sep_hist, 'Color', cmap(j,:), ...
             'DisplayName', sprintf('\\theta_2 = %.2f\\pi', th_vals(j)/pi), ...
             'HandleVisibility', vis);
    end
end
xlabel('{\it t}');
ylabel('separation');
legend('Location', 'best', 'Color', 'none');
ylim([0 2*R0])

% --- functions ---------------------------------------------------------------
function [t_hist, sep_hist] = simulate(x, th, beta, lam, uS, dt, tf)
    Nt       = round(tf / dt);
    t_hist   = zeros(Nt, 1);
    sep_hist = zeros(Nt, 1);
    for nt = 1:Nt
        [k1u, k1w] = rhs(x,             th,             beta, lam, uS);
        [k2u, k2w] = rhs(x + dt/2*k1u,  th + dt/2*k1w,  beta, lam, uS);
        [k3u, k3w] = rhs(x + dt/2*k2u,  th + dt/2*k2w,  beta, lam, uS);
        [k4u, k4w] = rhs(x + dt  *k3u,  th + dt  *k3w,  beta, lam, uS);
        x  = x  + dt/6 * (k1u + 2*k2u + 2*k3u + k4u);
        th = th + dt/6 * (k1w + 2*k2w + 2*k3w + k4w);
        t_hist(nt)   = nt * dt;
        sep_hist(nt) = norm(x(1,:) - x(2,:));
        if sep_hist(nt) < 5e-2; break; end
    end
    t_hist   = t_hist(1:nt);
    sep_hist = sep_hist(1:nt);
end


function u = uinf_pair(xi, xj, ~, pj)
    rvec = xi - xj;           % 1×2 row
    r    = norm(rvec);
    rhat = rvec / r;           % 1×2 row
    M    = (eye(2) - (rhat'*rhat)) / (4*pi*r);
    u    = M * pj';
end

function w = winf_pair(xi, xj, pi_, pj)
    rvec  = xi - xj;
    r     = norm(rvec);
    rhat  = rvec / r;
    rdotp = rhat * pj';       % scalar: r̂·p
    M     = ((eye(2) - 3*rhat'*rhat) * rdotp + rhat'*pj + pj'*rhat) / (4*pi*r^2);
    qi    = [-pi_(2), pi_(1)];
    w     = -(qi * M * pi_');
end

function [u, w] = rhs(x, th, beta, lam, uS)
    N = size(x, 1);
    p = [cos(th), sin(th)];
    u = uS * p;
    w = zeros(N, 1);
    for ni = 1:N
        for nj = 1:N
            if ni ~= nj
                u(ni,:) = u(ni,:) + beta      * uinf_pair(x(ni,:), x(nj,:), p(ni,:), p(nj,:))';
                w(ni)   = w(ni)   + beta * lam * winf_pair(x(ni,:), x(nj,:), p(ni,:), p(nj,:));
            end
        end
    end
end
