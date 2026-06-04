
close all
fig = journal_figure([6.5 2.25], 2, 'tex');

colormap = cmocean('tempo', 3);
c1 = colormap(3,:);
c2 = colormap(2,:);

% ── (a) Monopole coefficient vs separation ────────────────────────────────────

data = load('sigma_vs_separation.dat');
d        = data(:,1);
I1_th0   = data(:,2);
I1_thpi2 = data(:,4);
d_fine = linspace(d(1), d(end), 300)';

% Fit power-law exponent on second half of data
idx  = floor(numel(d)/2)+1 : numel(d);
p    = polyfit(log(d(idx)), log(abs(I1_th0(idx))), 1);
alpha    = p(1);
fit0 = exp(polyval(p, log(d_fine)));
fprintf('abs(I1_th0) ~ d^%.4f\n', alpha);

sp1 = subplot(1, 2, 1);
plot(d, abs(I1_th0), 'Color', c1, 'DisplayName', '{\it \theta} = 0'); hold on
plot(d, abs(I1_thpi2), 'Color', c2, 'DisplayName', '{\it \theta} = {\it \pi}/2');
plot(d_fine, fit0, '--', 'Color', [0.2 0.2 0.2], ...
     'DisplayName', sprintf('\\propto ({\\it \\Delta x}_{cm})^{%.2f}', alpha));
xlabel('separation, {\it \Delta x}_{cm}')
ylabel('monopole, {\it m}_0')
legend('Location', 'northeast', 'Color', 'none', 'FontSize', sp1.FontSize)
xlim([0 100])
ylim([0 0.1])

% ── (b) Multipole expansion error vs r ───────────────────────────────────────

err = load('multipole_error_2p.dat');
r         = err(:,1);
err_th0   = abs(err(:,2));
err_thpi2 = abs(err(:,3));

sp2 = subplot(1, 2, 2);
loglog(r, err_th0, '-',  'Color', c1, 'DisplayName', '{\it \theta} = 0'); hold on
loglog(r, err_thpi2, '-', 'Color', c2, 'DisplayName', '{\it \theta} = {\it \pi}/2');
loglog(r, 5 * err_th0(1) * (r(1) ./ r).^1, 'k--', 'DisplayName', sprintf('\\propto {\\it r}^{-1}'));
loglog(r, (1/5) * err_thpi2(1) * (r(1) ./ r).^2, 'k:',   'DisplayName', sprintf('\\propto {\\it r}^{-2}'));
xlabel('{\it r}')
ylabel('|{\bf\it u} - {\bf\it u}_{mp}|')
legend('Location', 'southwest', 'Color', 'none', 'FontSize', sp2.FontSize)
% xlim([1e2 1e4])
ylim([1e-12 1e-2])
yticks(10.^(-12:2:2))

% ── Layout ────────────────────────────────────────────────────────────────────

sp1.Units = 'inches';
sp2.Units = 'inches';

sp1.Position(1) = 1;
sp1.Position(2) = 0.9;
sp1.Position(3) = 0.4 * fig.PaperSize(1);
sp1.Position(4) = 0.75 * fig.PaperSize(2);

sp2.Position(1) = sum(sp1.Position([1 3])) + 1.25;
sp2.Position(2:4) = sp1.Position(2:4);

a = annotation('textbox', 'string', '{\it (a)}', 'FontName', 'Times', 'EdgeColor', 'none', 'FontSize', sp1.FontSize, 'Units', 'inches');
a.Position(1) = sp1.Position(1)-1;
a.Position(2) = sum(sp1.Position([2 4])) - 0.5*a.Position(4);

b = annotation('textbox', 'string', '{\it (b)}', 'FontName', 'Times', 'EdgeColor', 'none', 'FontSize', sp2.FontSize, 'Units', 'inches');
b.Position(1) = sp2.Position(1)-1.1;
b.Position(2) = sum(sp2.Position([2 4])) - 0.5*b.Position(4);
