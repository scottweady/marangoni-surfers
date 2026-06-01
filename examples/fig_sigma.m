
close all
fig = journal_figure([6.5 2], 2, 'tex');

c1 = [0.2 0.2 0.8];
c2 = [0.8 0.2 0.8];

% ── (a) Monopole coefficient vs separation ────────────────────────────────────

data = load('sigma_vs_separation.dat');
d        = data(:,1);
I1_th0   = data(:,2);
I1_thpi2 = data(:,4);
d_fine = linspace(d(1), d(end), 300)';

sp1 = subplot(1, 2, 1);
plot(d, abs(I1_th0), 'Color', c1, 'DisplayName', '{\it \theta} = 0'); hold on
plot(d, abs(I1_thpi2), 'Color', c2, 'DisplayName', '{\it \theta} = {\it \pi}/2');
plot(d_fine, abs(I1_th0(1)) ./ (d_fine/d_fine(1)).^(3/2), '--', 'Color', [0.2 0.2 0.2], ...
     'DisplayName', '\propto {\it d}^{-3/2}');
xlabel('separation, {\it d}')
ylabel('monopole, {\it m}_0')
legend('Location', 'northeast', 'Color', 'none', 'FontSize', sp1.FontSize)
% ylim([-0.05 0.25])
ylim auto

% ── (b) Multipole expansion error vs r ───────────────────────────────────────

err = load('multipole_error_2p.dat');
r         = err(:,1);
err_th0   = err(:,2);
err_thpi2 = err(:,3);

sp2 = subplot(1, 2, 2);
loglog(r, err_th0, 'o-',  'Color', c1, 'DisplayName', '{\it \theta} = 0'); hold on
loglog(r, err_thpi2, 'o-', 'Color', c2, 'DisplayName', '{\it \theta} = {\it \pi}/2');
% loglog(r, 5 * err_th0(1) * (r(1) ./ r).^1, 'k--',   'HandleVisibility', 'off');
loglog(r, (1/5) * err_thpi2(1) * (r(1) ./ r).^2, 'k:',   'HandleVisibility', 'off');
xlabel('{\it r}')
ylabel('|{\bf u}_{num} - {\bf u}_{mp}|')
legend('Location', 'southwest', 'Color', 'none', 'FontSize', sp2.FontSize)
% xlim([1e2 1e4])
% ylim([1e-8 1e-2])

% ── Layout ────────────────────────────────────────────────────────────────────

sp1.Units = 'inches';
sp2.Units = 'inches';

sp1.Position(1) = 1;
sp1.Position(2) = 0.75;
sp1.Position(3) = 0.4 * fig.PaperSize(1);
sp1.Position(4) = 0.75 * fig.PaperSize(2);

sp2.Position(1) = sum(sp1.Position([1 3])) + 1.25;
sp2.Position(2:4) = sp1.Position(2:4);

a = annotation('textbox', 'string', '{\it (a)}', 'FontName', 'Times', 'EdgeColor', 'none', 'FontSize', 18, 'Units', 'inches');
a.Position(1) = sp1.Position(1);
a.Position(2) = sum(sp1.Position([2 4])) - a.Position(4);

b = annotation('textbox', 'string', '{\it (b)}', 'FontName', 'Times', 'EdgeColor', 'none', 'FontSize', 18, 'Units', 'inches');
b.Position(1) = sp2.Position(1);
b.Position(2) = sum(sp2.Position([2 4])) - b.Position(4);
