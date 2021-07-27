clear;

deg = 1.11;

figure('Position', [1 1 .6 * 1920 .35 * 1080]); %[1 1 .6 * 1920 .5 * 1080]

strDisp{1} = '$n=-10^{11} \mathrm{cm}^{-2}$';
strDisp{2} = '$n=-2\times 10^{11} \mathrm{cm}^{-2}$';
strDisp{3} = '$n=-3\times 10^{11} \mathrm{cm}^{-2}$';
strDisp{4} = '$n=-4\times 10^{11} \mathrm{cm}^{-2}$';

xtickvec = [0 50 100];
ytickvec = 0.1:0.1:0.5;

xlimvec = [0 160];
ylimvec = [0 0.5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Theory
ntheory = [10 20 30 40];
Ttheory = linspace(.1, 160, 30);
teV = 3.4;

%Theory e-ph
load(['parameph_deg' num2str(deg) '_e.mat']);
paramcell = num2cell(param);
paramdeltacell = num2cell(paramdelta);
[betaAtildeeVe, vbare] = paramcell{:};
[deltabetaAtildeeVe, deltavbare] = paramdeltacell{:};

load(['parameph_deg' num2str(deg) '_h.mat']);
paramcell = num2cell(param);
paramdeltacell = num2cell(paramdelta);
[betaAtildeeVh, vbarh] = paramcell{:};
[deltabetaAtildeeVh, deltavbarh] = paramdeltacell{:};

clr = jet(numel(ntheory));

%subplot(1, 3, 1);
subplot('Position', [0.08 .2 0.29 .75]);
hold on;
box on;

for j = 1:4 %coz Dirac model is weakly density dependent anyway
    n = ntheory(j);

    load(['rhoeph_vsT_deg' num2str(deg) '_n' num2str(n) '_Dirac_e.mat'], 'Ttheory', 'n', 'deg', 'betaAtildeeVe', 'vbare', 'rhoeph_Dirac_e');
    load(['rhoeph_vsT_deg' num2str(deg) '_n' num2str(n) '_Dirac_h.mat'], 'Ttheory', 'n', 'deg', 'betaAtildeeVh', 'vbarh', 'rhoeph_Dirac_h');
    load(['rhoeph_vsT_deg' num2str(deg) '_n' num2str(n) '_CN_e.mat'], 'Ttheory', 'n', 'deg', 'betaAtildeeVe', 'vbare', 'rhoeph_CN_e');
    load(['rhoeph_vsT_deg' num2str(deg) '_n' num2str(n) '_CN_h.mat'], 'Ttheory', 'n', 'deg', 'betaAtildeeVh', 'vbarh', 'rhoeph_CN_h');        
    
    index = rhoeph_Dirac_h < 0.417;
    p(j) = plot(Ttheory(index), rhoeph_Dirac_h(index), '--', 'Color', clr(j, :), 'LineWidth', 3, 'DisplayName', sprintf('Electron-phonon\nDirac model')); %['$n=' sprintf('%.0f', abs(n)) '\times 10^{10}\mathrm{cm}^{-2}$']

    plot(Ttheory, rhoeph_CN_h, '-', 'Color', clr(j, :), 'LineWidth', 3, 'DisplayName', sprintf('Electron-phonon\nCastro-Neto model'));
end


xticks(xtickvec);
xticklabels(num2cell(abs(xtickvec)));

xlabel('$T$~(K)', 'FontSize', 30, 'Interpreter', 'latex');
ylabel('$\rho~(h/e^2)$', 'FontSize', 30, 'Interpreter', 'latex');

annotation('textbox', [0.1578 0.8439 0.1329 0.0761], 'String', 'Phonon theory', 'Interpreter', 'latex', 'FontSize', 22, 'LineWidth', 2); %, 'LineStyle', 'none'

strparamh{1} = ['$\tilde{\beta}_A^{(h)}=' sprintf('%.0f', betaAtildeeVh) ' \pm ' sprintf('%.0f', deltabetaAtildeeVh) '$~eV'];
strparamh{2} = ['$v_F^{(h)}/v_0=' sprintf('%.2f', vbarh) ' \pm ' sprintf('%.2f', deltavbarh) '$'];
annotation('textbox', [0.22 0.25  .1 .1], 'String', strparamh, 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 17);

set(gca, 'LineWidth', 3, 'FontSize', 25, 'FontWeight', 'bold');

xlim(xlimvec);
ylim(ylimvec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexp = -ntheory;
clr = jet(length(nexp));

%subplot(1, 3, 2);
subplot('Position', [0.37 .2 0.29 .75]); %[0.37 .13 0.29 .8]
hold on;
box on;

Tmin = 0;
    
for j = 1:length(nexp)
    n = nexp(j);
    disp(['n=' num2str(n)]);

    load(['rho_experiment_vsT_deg' num2str(deg) '_n' num2str(n) '.mat'], 'n', 'deg', 'Texp', 'rhoexp');
    p(j) = plot(Texp, rhoexp, 's', 'Color', clr(j, :), 'MarkerSize', 3, 'LineWidth', 2, 'DisplayName', strDisp{j});
end

rectangle('Position', [1.6 0.006 20 0.1], 'LineWidth', 2, 'LineStyle', '--');
annotation('textarrow', [0.4345 0.4089], [0.4280 0.3571], 'String', 'see Fig.~7a', 'FontSize', 15, 'Interpreter', 'latex');

xlabel('$T$~(K)', 'FontSize', 30, 'Interpreter', 'latex');

xticks(xtickvec);
xticklabels(num2cell(abs(xtickvec)));

yticks(ytickvec);
yticklabels({});

l = legend(p);
legend boxoff;
set(l, 'Interpreter', 'latex', 'Location', 'SouthEast', 'FontSize', 17);

annotation('textbox', [.38 .73 .1 .1], 'String', ['$\theta=' num2str(deg) '^\circ$'], 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 22);
annotation('textbox', [0.4248 0.8307 0.1889 0.0868], 'String', 'Experiment (Ref.~13)', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 22);

set(gca, 'LineWidth', 3, 'FontSize', 25, 'FontWeight', 'bold');

xlim(xlimvec);
ylim(ylimvec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Theory Planckian
load(['paramPlanckian_deg' num2str(deg) '_e.mat']);
paramcell = num2cell(param);
paramdeltacell = num2cell(paramdelta);
[Ce, vbare] = paramcell{:};
[deltaCe, deltavbare] = paramdeltacell{:};

load(['paramPlanckian_deg' num2str(deg) '_h.mat']);
paramcell = num2cell(param);
paramdeltacell = num2cell(paramdelta);
[Ch, vbarh] = paramcell{:};
[deltaCh, deltavbarh] = paramdeltacell{:};

clr = jet(length(ntheory));

%subplot(1, 3, 3);
subplot('Position', [0.66 .2 0.29 .75]);
hold on;
box on;

for j = 1:length(ntheory)
    n = ntheory(j);
   
    load(['rho_Planckian_Dirac_deg' num2str(deg) '_n' sprintf('%.2f', n) '_e.mat'], 'Ttheory', 'n', 'deg', 'teV', 'Ce', 'vbare', 'rhoDirac_e');
    load(['rho_Planckian_Dirac_deg' num2str(deg) '_n' sprintf('%.2f', n) '_h.mat'], 'Ttheory', 'n', 'deg', 'teV', 'Ch', 'vbarh', 'rhoDirac_h');

    load(['rho_Planckian_CN_deg' num2str(deg) '_n' sprintf('%.2f', n) '_e.mat'], 'Ttheory', 'n', 'deg', 'teV', 'Ce', 'vbare', 'rhoe');
    load(['rho_Planckian_CN_deg' num2str(deg) '_n' sprintf('%.2f', n) '_h.mat'], 'Ttheory', 'n', 'deg', 'teV', 'Ch', 'vbarh', 'rhoh');

    p(j) = plot(Ttheory, rhoh, '-', 'Color', clr(j, :), 'LineWidth', 3, 'DisplayName', strDisp{j}); %, , 'DisplayName', strDisp{j}
    
    plot(Ttheory, rhoDirac_h, '--', 'Color', clr(j, :), 'LineWidth', 3);
end

xticks(xtickvec);
xticklabels(num2cell(abs(xtickvec)));

yticks(ytickvec);
yticklabels({});

xlabel('$T$~(K)', 'FontSize', 30, 'Interpreter', 'latex');

annotation('textbox', [0.7313 0.8439 0.1515 0.0761], 'String', 'Planckian theory', 'Interpreter', 'latex', 'FontSize', 22, 'LineWidth', 2);

strparamh{1} = ['$C_h=' sprintf('%.1f', Ch) ' \pm ' sprintf('%.1f', deltaCh) '$'];
strparamh{2} = ['$v_F^{(h)}/v_0=' sprintf('%.2f', vbarh) ' \pm ' sprintf('%.2f', deltavbarh) '$'];
annotation('textbox', [0.8 0.22 .1 .1], 'String', strparamh, 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 17);

x0 = .09;
dx = 0.29;
y = .81;
annotation('textbox', [0 y .1 .1], 'String', '(b)', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 25);

set(gca, 'LineWidth', 3, 'FontSize', 25, 'FontWeight', 'bold');

xlim(xlimvec);
ylim(ylimvec);

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
print('-dpdf', ['rho_vsT_beyondDirac_deg' num2str(deg) '.pdf']);

% Restore the previous settings
set(gcf,'PaperType',prePaperType);
set(gcf,'PaperUnits',prePaperUnits);
set(gcf,'Units',preUnits);
set(gcf,'PaperPosition',prePaperPosition);
set(gcf,'PaperSize',prePaperSize);