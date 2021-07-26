clear;

deg = 1.24;

%Theory
degstr = num2str(deg);
degstr = strrep(degstr, '.', '_');

figure();
hold on;
box on;

j = 1;
n = -0.02;

load(['rho_tot_emt_deg' num2str(deg) '_n' sprintf('%.2f', n) '_CN.mat'], 'Ttheory', 'n', 'deg', 'sigmas', 'nimp', 'nrms', 'alpha0', 'betaAeV', 'rhotot');
        
if sprintf('%.0f', abs(n)) == '0'
    txtunit = '';
else
    txtunit = '\times 10^{10}\mathrm{cm}^{-2}';
end    
    
p(j) = plot(Ttheory, rhotot, '-', 'Color', [0 0 0], 'LineWidth', 3);

T1 = 3.5;
T2 = 16;
T3 = 65;
ybound = [0.0 0.145];
plot(T1 * [1 1], ybound, 'k--', 'LineWidth', 2);
plot(T2 * [1 1], ybound, 'k--', 'LineWidth', 2);
plot(T3 * [1 1], ybound, 'k--', 'LineWidth', 2);

xlabel('$T~(K)$', 'FontSize', 30, 'Interpreter', 'latex');
ylabel('$\rho~(h/e^2)$', 'FontSize', 30, 'Interpreter', 'latex');

annotation('textbox', [0.28 0.82 0.2342 0.0952], 'String', '$T \gg (T_F,T_\mathrm{BG}) $', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 22);
annotation('textbox', [0.28 0.75 0.2342 0.0952], 'String', 'phonons', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 22);

annotation('textbox', [0.5629 0.82 0.2342 0.0952], 'String', '$T \sim \varepsilon_\mathrm{VHS} $', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 22);
annotation('textbox', [0.5629 0.75 0.2342 0.0952], 'String', 'phonons', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 22);

x0 = 1;
y0 = 0.07;
w = 7;
h = 0.017;
rectangle('Position',[x0 y0 w h]);
plot([x0 41], [y0 0.021], 'k');
plot([(x0 + w) 122], [(y0 + h) 0.071], 'k');

xlim([0 125]);
ylim(ybound);

set(gca, 'LineWidth', 3, 'FontSize', 25, 'FontWeight', 'bold');

%Inset%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes('Position',[.41 .28 .48 .26])
box on
hold on

plot(Ttheory, rhotot, '-', 'Color', [0 0 0], 'LineWidth', 3);

xlabel('$T~(K)$', 'FontSize', 30, 'Interpreter', 'latex');
ylabel('$\rho~(h/e^2)$', 'FontSize', 30, 'Interpreter', 'latex');

ybound = [y0 (y0 + h)];
plot(T1 * [1 1], ybound, 'k--', 'LineWidth', 2);

annotation('textbox', [0.4157 0.3057 0.2342 0.0952], 'String', '$T \ll T_F$', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 18);
annotation('textbox', [0.4139 0.2595 0.2342 0.0952], 'String', 'impurity', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 18);

% annotation('textbox', [0.59 0.3057 0.2342 0.0952], 'String', '$T > T_F$', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 18);
annotation('textbox', [0.59 0.3057 0.2342 0.0952], 'String', '$T$', 'LineStyle', 'none', 'Interpreter', 'latex', 'FontSize', 18);
annotation('textbox', [0.62 0.328 0.2342 0.0952], 'String', char(8819), 'LineStyle', 'none', 'FontSize', 29);
annotation('textbox', [0.66 0.3057 0.2342 0.0952], 'String', '$T_F$', 'LineStyle', 'none', 'Interpreter', 'latex', 'FontSize', 18);
annotation('textbox', [0.59 0.2595 0.2342 0.0952], 'String', 'impurity', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 18);

xlim([x0 (x0 + w)]);
ylim(ybound);

set(gca, 'LineWidth', 2, 'FontSize', 15, 'FontWeight', 'bold');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print('-dpdf', 'regime.pdf');