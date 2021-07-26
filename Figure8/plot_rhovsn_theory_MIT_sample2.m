clear;

%Units
%T = K
%n = 10^10 cm^-2

deg = 1.07;

%Experiment
filename = ['resistivity_experiment_deg' num2str(deg) '.mat'];
load(filename, 'n_exp', 'rho_exp');

%Phonon
filename = ['resistivity_phonon_deg' num2str(deg) '.mat'];
load(filename, 'T', 'n', 'deg', 'teV', 'betaAtilde', 'vbar', 'rho_phonon');

%Planckian
filename = ['resistivity_Planckian_deg' num2str(deg) '.mat'];
load(filename, 'T', 'n', 'deg', 'teV', 'C', 'vbar', 'rho_Planckian');

%Planckian EMT
nrms = 5;
filename = ['resistivity_PlanckianEMT_deg' num2str(deg) '.mat'];
load(filename, 'T', 'n', 'deg', 'teV', 'C', 'vbar', 'rho_PlanckianEMT');

figure();
hold on;
box on;

p(1) = plot(n_exp, (rho_exp - rho_exp(1)) / rho_exp(1), 'ks', 'MarkerSize', 7, 'MarkerFaceColor', 'k', 'DisplayName', 'Experiment (MIT, Ref.~60)');
p(2) = plot(n, (rho_phonon - rho_phonon(1)) /rho_phonon(1) , 'b', 'Linewidth', 3, 'DisplayName', 'Phonon'); 
p(3) = plot(n, (rho_Planckian - rho_Planckian(1)) / rho_Planckian(1), 'r', 'Linewidth', 3, 'DisplayName', 'Planckian'); 
p(4) = plot(n, (rho_PlanckianEMT - rho_PlanckianEMT(1)) / rho_PlanckianEMT(1), 'g', 'Linewidth', 3, 'DisplayName', 'Planckian EMT');

xlabel('$n~(10^{10}\mathrm{cm}^{-2})$', 'FontSize', 30, 'Interpreter', 'latex');
ylabel('$[\rho-\rho(0)]/\rho(0)$', 'FontSize', 30, 'Interpreter', 'latex');

l = legend(p);
legend boxoff;
set(l, 'Interpreter', 'latex', 'FontSize', 16, 'Position', [.65 .6 .1 .1]); %, 'Position', [.385 .48 .1 .1]

annotation('textbox', [.18 .29 .1 .1], 'String', {['$\theta=' num2str(deg) '^\circ$'],['$T = ' num2str(T) '~\mathrm{K}$'], ['$n_\mathrm{rms}=' num2str(nrms) '\times 10^{10}\mathrm{cm}^{-2}$']}, 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 25);
annotation('textbox', [.72 .15 .1 .1], 'String', 'Electron', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 25);
% annotation('textbox', [.49 .17 .1 .1], 'String', ['$n_\mathrm{rms}=' num2str(nrms) '\times 10^{10}\mathrm{cm}^{-2}$'], 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 25);

set(gca, 'LineWidth', 3, 'FontSize', 25, 'FontWeight', 'bold');

xlim([0 40]);
ylim([-0.5 0.12]);

%Save
% Backup previous settings
prePaperType = get(gcf,'PaperType');
prePaperUnits = get(gcf,'PaperUnits');
preUnits = get(gcf,'Units');
prePaperPosition = get(gcf,'PaperPosition');
prePaperSize = get(gcf,'PaperSize');

% Make changing paper type possible
set(gcf,'PaperType','<custom>');

% Set units to all be the same
set(gcf,'PaperUnits','inches');
set(gcf,'Units','inches');

% Set the page size and position to match the figure's dimensions
%paperPosition = get(gcf,'PaperPosition');
position = get(gcf,'Position');
set(gcf,'PaperPosition',[0,0,position(3:4)]);
set(gcf,'PaperSize',position(3:4));

% Save the pdf
print('-dpdf', ['rhovsn_deg' num2str(deg) '.pdf']);

% Restore the previous settings
set(gcf,'PaperType',prePaperType);
set(gcf,'PaperUnits',prePaperUnits);
set(gcf,'Units',preUnits);
set(gcf,'PaperPosition',prePaperPosition);
set(gcf,'PaperSize',prePaperSize);    