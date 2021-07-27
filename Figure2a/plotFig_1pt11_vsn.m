clear;

deg = 1.11;

figure('Position', [1 1 .6 * 1920 .35 * 1080]); %[1 1 .6 * 1920 .5 * 1080]

strDisp = @(T1) ['$T=' sprintf('%.0f', T1) '$~K'];

xtickvec = -20:20:20;
ytickvec = 0.1:0.1:0.5;

xlimvec = [-1 1] * 40;
ylimvec = [0 0.5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Theory
Ttheory = [50 80 100 150];
nmax = 40; %nVHSfunCN(deg)
ntheory = linspace(0.1, nmax, 30);
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

clr = jet(numel(Ttheory));

%subplot(1, 3, 1);
subplot('Position', [0.08 .2 0.29 .75]);
hold on;
box on;

for j = 1:numel(Ttheory) %coz Dirac model is weakly density dependent anyway
    T = Ttheory(j);    

    load(['rhoeph_vsn_deg' num2str(deg) '_T' num2str(T) '_Dirac_e.mat'], 'ntheory', 'T', 'deg', 'betaAtildeeVe', 'vbare', 'rhoeph_Dirac_e');
    load(['rhoeph_vsn_deg' num2str(deg) '_T' num2str(T) '_Dirac_h.mat'], 'ntheory', 'T', 'deg', 'betaAtildeeVh', 'vbarh', 'rhoeph_Dirac_h');
    load(['rhoeph_vsn_deg' num2str(deg) '_T' num2str(T) '_CN_e.mat'], 'ntheory', 'T', 'deg', 'betaAtildeeVe', 'vbare', 'rhoeph_CN_e');
    load(['rhoeph_vsn_deg' num2str(deg) '_T' num2str(T) '_CN_h.mat'], 'ntheory', 'T', 'deg', 'betaAtildeeVh', 'vbarh', 'rhoeph_CN_h');
        
    p(j) = plot([fliplr(-ntheory) NaN ntheory], [fliplr(rhoeph_Dirac_h) NaN rhoeph_Dirac_e], '--', 'Color', clr(j, :), 'LineWidth', 3, 'DisplayName', sprintf('Electron-phonon\nDirac model')); %['$n=' sprintf('%.0f', abs(n)) '\times 10^{10}\mathrm{cm}^{-2}$']
    plot([fliplr(-ntheory) NaN ntheory], [fliplr(rhoeph_CN_h) NaN rhoeph_CN_e], '-', 'Color', clr(j, :), 'LineWidth', 3, 'DisplayName', sprintf('Electron-phonon\nCastro-Neto model'));
end

plot([0 0], 0.86 * ylimvec, 'k--', 'LineWidth', 3);

xticks(xtickvec);
xticklabels(num2cell(xtickvec));

xlabel('$n~(10^{10}\mathrm{cm}^{-2})$', 'FontSize', 30, 'Interpreter', 'latex');
ylabel('$\rho~(h/e^2)$', 'FontSize', 30, 'Interpreter', 'latex');

annotation('textbox', [.085 .18 .1 .1], 'String', 'Holes', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 20);
annotation('textbox', [.292 .18 .1 .1], 'String', 'Electron', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 20);
a = annotation('textbox', [0.1578 0.8439 0.1330 0.0861], 'String', 'Phonon theory','Interpreter', 'latex', 'FontSize', 22, 'LineWidth', 2);

strparame{1} = ['$\tilde{\beta}_A^{(e)}=' sprintf('%.0f', betaAtildeeVe) ' \pm ' sprintf('%.0f', deltabetaAtildeeVe) '$~eV'];
strparame{2} = ['$v_F^{(e)}/v_0=' sprintf('%.2f', vbare) ' \pm ' sprintf('%.2f', deltavbare) '$'];
annotation('textbox', [.24 .34 .1 .1], 'String', strparame, 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 15);
strparamh{1} = ['$\tilde{\beta}_A^{(h)}=' sprintf('%.0f', betaAtildeeVh) ' \pm ' sprintf('%.0f', deltabetaAtildeeVh) '$~eV'];
strparamh{2} = ['$v_F^{(h)}/v_0=' sprintf('%.2f', vbarh) ' \pm ' sprintf('%.2f', deltavbarh) '$'];
annotation('textbox', [.082 .34 .1 .1], 'String', strparamh, 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 15);

set(gca, 'LineWidth', 3, 'FontSize', 25, 'FontWeight', 'bold');

xlim(xlimvec);
ylim(ylimvec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subplot(1, 3, 2);
subplot('Position', [0.37 .2 0.29 .75]);
hold on;
box on;

Texp = Ttheory;
nmin = 0;
    
for j = 1:length(Texp)
    T = Texp(j);    
    load(['rho_experiment_vsn_deg' num2str(deg) '_T' num2str(T) '.mat'], 'T', 'deg', 'nexp', 'rhoexp');
    
    p(j) = plot(nexp, rhoexp, 's', 'Color', clr(j, :), 'MarkerSize', 3, 'LineWidth', 2, 'DisplayName', strDisp(T)); %'MarkerFaceColor', clr(j, :)
end

plot([0 0], 0.86 * ylimvec, 'k--', 'LineWidth', 3);

xlabel('$n~(10^{10}\mathrm{cm}^{-2})$', 'FontSize', 30, 'Interpreter', 'latex');
%ylabel('$\rho~(h/e^2)$', 'FontSize', 30, 'Interpreter', 'latex');

xticks(xtickvec);
xticklabels(num2cell(xtickvec));

yticks(ytickvec);
yticklabels({});

l = legend(p);
legend boxoff;
set(l, 'Interpreter', 'latex', 'Position', [.38 .38 .1 .1], 'FontSize', 20);

annotation('textbox', [.375 .765 .1 .1], 'String', ['$\theta=' num2str(deg) '^\circ$'], 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 22);
annotation('textbox', [.375 .18 .1 .1], 'String', 'Holes', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 20);
annotation('textbox', [.582 .18 .1 .1], 'String', 'Electron', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 20);
annotation('textbox', [.43 .83 .1 .1], 'String', 'Experiment (Ref.~13)', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 22);

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

clr = jet(length(Ttheory));

%subplot(1, 3, 3);
subplot('Position', [0.66 .2 0.29 .75]);
hold on;
box on;

for j = 1:length(Ttheory)
    T = Ttheory(j);
   
    load(['rho_Planckian_Dirac_deg' num2str(deg) '_T' num2str(T) '_e.mat'], 'T', 'ntheory', 'deg', 'teV', 'Ce', 'vbare', 'rhoDirac_e');
    load(['rho_Planckian_Dirac_deg' num2str(deg) '_T' num2str(T) '_h.mat'], 'T', 'ntheory', 'deg', 'teV', 'Ch', 'vbarh', 'rhoDirac_h');

    load(['rho_Planckian_CN_deg' num2str(deg) '_T' num2str(T) '_e.mat'], 'T', 'ntheory', 'deg', 'teV', 'Ce', 'vbare', 'rhoe');
    load(['rho_Planckian_CN_deg' num2str(deg) '_T' num2str(T) '_h.mat'], 'T', 'ntheory', 'deg', 'teV', 'Ch', 'vbarh', 'rhoh');
    
    p(j) = plot([fliplr(-ntheory) ntheory], [fliplr(rhoh) rhoe], '-', 'Color', clr(j, :), 'LineWidth', 3, 'DisplayName', strDisp(T)); %, , 'DisplayName', strDisp{j}
    
    plot([fliplr(-ntheory) ntheory], [fliplr(rhoDirac_h) rhoDirac_e], '--', 'Color', clr(j, :), 'LineWidth', 3);  
end

plot([0 0], 0.86 * ylimvec, 'k--', 'LineWidth', 3); %'Color', [.4 .4 .4]

xticks(xtickvec);
xticklabels(num2cell(xtickvec));

yticks(ytickvec);
yticklabels({});

xlabel('$n~(10^{10}\mathrm{cm}^{-2})$', 'FontSize', 30, 'Interpreter', 'latex');

annotation('textbox', [.665 .18 .1 .1], 'String', 'Holes', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 20);
annotation('textbox', [.872 .18 .1 .1], 'String', 'Electron', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 20);
annotation('textbox', [0.7296 0.8439 0.1524 0.0847], 'String', 'Planckian theory', 'Interpreter', 'latex', 'FontSize', 22, 'LineWidth', 2);

strparame{1} = ['$C_e=' sprintf('%.1f', Ce) ' \pm ' sprintf('%.1f', deltaCe) '$'];
strparame{2} = ['$v_F^{(e)}/v_0=' sprintf('%.1f', vbare) ' \pm ' sprintf('%.1f', deltavbare) '$'];
annotation('textbox', [.835 .3 .1 .1], 'String', strparame, 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 15);
strparamh{1} = ['$C_h=' sprintf('%.1f', Ch) ' \pm ' sprintf('%.1f', deltaCh) '$'];
strparamh{2} = ['$v_F^{(h)}/v_0=' sprintf('%.2f', vbarh) ' \pm ' sprintf('%.2f', deltavbarh) '$'];
annotation('textbox', [.66 .3 .1 .1], 'String', strparamh, 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 15);

x0 = .09;
dx = 0.29;
y = .83;
annotation('textbox', [0 y .1 .1], 'String', '(a)', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 25);

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
print('-dpdf', ['rho_vsn_beyondDirac_deg' num2str(deg) '.pdf']);

% Restore the previous settings
set(gcf,'PaperType',prePaperType);
set(gcf,'PaperUnits',prePaperUnits);
set(gcf,'Units',preUnits);
set(gcf,'PaperPosition',prePaperPosition);
set(gcf,'PaperSize',prePaperSize);    