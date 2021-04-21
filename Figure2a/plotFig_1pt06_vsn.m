%function plotFig_1pt06()

clear;

savefile = 1;
loadfile = 1;

deg = 1.06;

figure('Position', [1 1 .6 * 1920 .35 * 1080]); %[1 1 .6 * 1920 .5 * 1080]

strDisp = @(T1) ['$T=' sprintf('%.0f', T1) '$~K'];

xtickvec = -20:20:20;
ytickvec = 0.2:0.2:0.6;

xlimvec = [-1 1] * 40;
ylimvec = [0 0.6];

xshift = 0.02;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Theory e-ph
Ttheory = [50 80 100 150];
nmax = 40;
ntheory = linspace(.1, nmax, 30);
teV = 3.4;

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
subplot('Position', [xshift + 0.08 .2 0.29 .75]);
hold on;
box on;

for j = 1:numel(Ttheory) %coz Dirac model is weakly density dependent anyway
    T = Ttheory(j);
        
    if loadfile
        load(['rhoeph_vsn_deg' num2str(deg) '_T' num2str(T) '_Dirac_e.mat'], 'ntheory', 'T', 'deg', 'betaAtildeeVe', 'vbare', 'rhoeph_Dirac_e');
        load(['rhoeph_vsn_deg' num2str(deg) '_T' num2str(T) '_Dirac_h.mat'], 'ntheory', 'T', 'deg', 'betaAtildeeVh', 'vbarh', 'rhoeph_Dirac_h');
        load(['rhoeph_vsn_deg' num2str(deg) '_T' num2str(T) '_CN_e.mat'], 'ntheory', 'T', 'deg', 'betaAtildeeVe', 'vbare', 'rhoeph_CN_e');
        load(['rhoeph_vsn_deg' num2str(deg) '_T' num2str(T) '_CN_h.mat'], 'ntheory', 'T', 'deg', 'betaAtildeeVh', 'vbarh', 'rhoeph_CN_h');
    else
        rhoeph_Dirac_e = 1 ./ sigmaephfun(T, ntheory, deg, teV, betaAtildeeVe, vbare);
        rhoeph_Dirac_h = 1 ./ sigmaephfun(T, -ntheory, deg, teV, betaAtildeeVh, vbarh);
        
        rhoeph_CN_e = 1./ sigmaephCN_exact_avg_iso(T, ntheory, deg, [], teV, betaAtildeeVe, vbare);
        rhoeph_CN_h = 1./ sigmaephCN_exact_avg_iso(T, -ntheory, deg, [], teV, betaAtildeeVh, vbarh);
    end
    
    if savefile
        save(['rhoeph_vsn_deg' num2str(deg) '_T' num2str(T) '_Dirac_e.mat'], 'ntheory', 'T', 'deg', 'betaAtildeeVe', 'vbare', 'rhoeph_Dirac_e');
        save(['rhoeph_vsn_deg' num2str(deg) '_T' num2str(T) '_Dirac_h.mat'], 'ntheory', 'T', 'deg', 'betaAtildeeVh', 'vbarh', 'rhoeph_Dirac_h');
        save(['rhoeph_vsn_deg' num2str(deg) '_T' num2str(T) '_CN_e.mat'], 'ntheory', 'T', 'deg', 'betaAtildeeVe', 'vbare', 'rhoeph_CN_e');
        save(['rhoeph_vsn_deg' num2str(deg) '_T' num2str(T) '_CN_h.mat'], 'ntheory', 'T', 'deg', 'betaAtildeeVh', 'vbarh', 'rhoeph_CN_h');
    end
    
    p(j) = plot([fliplr(-ntheory) NaN ntheory], [fliplr(rhoeph_Dirac_h) NaN rhoeph_Dirac_e], '--', 'Color', clr(j, :), 'LineWidth', 3, 'DisplayName', sprintf('Electron-phonon\nDirac model')); %['$n=' sprintf('%.0f', abs(n)) '\times 10^{10}\mathrm{cm}^{-2}$']
    plot([fliplr(-ntheory) NaN ntheory], [fliplr(rhoeph_CN_h) NaN rhoeph_CN_e], '-', 'Color', clr(j, :), 'LineWidth', 3, 'DisplayName', sprintf('Electron-phonon\nCastro-Neto model'));
end

plot([0 0], 0.88 * ylimvec, 'k--', 'LineWidth', 3);

xticks(xtickvec);
xticklabels(num2cell(xtickvec));

xlabel('$n~(10^{10}\mathrm{cm}^{-2})$', 'FontSize', 30, 'Interpreter', 'latex');
ylabel('$\rho~(h/e^2)$', 'FontSize', 30, 'Interpreter', 'latex');

annotation('textbox', [xshift + .085 .18 .1 .1], 'String', 'Holes', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 20);
annotation('textbox', [xshift + .292 .18 .1 .1], 'String', 'Electron', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 20);
annotation('textbox', [0.1761 0.8545 0.1355 0.0755], 'String', 'Phonon theory','Interpreter', 'latex', 'FontSize', 22, 'LineWidth', 2);

strparame{1} = ['$\tilde{\beta}_A^{(e)}=' sprintf('%.0f', betaAtildeeVe) ' \pm ' sprintf('%.0f', deltabetaAtildeeVe) '$~eV'];
strparame{2} = ['$v_F^{(e)}/v_0=' sprintf('%.3f', vbare) ' \pm ' sprintf('%.3f', deltavbare) '$'];
annotation('textbox', [xshift + .23 .3 .1 .1], 'String', strparame, 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 15);
strparamh{1} = ['$\tilde{\beta}_A^{(h)}=' sprintf('%.0f', betaAtildeeVh) ' \pm ' sprintf('%.0f', deltabetaAtildeeVh) '$~eV'];
strparamh{2} = ['$v_F^{(h)}/v_0=' sprintf('%.2f', vbarh) ' \pm ' sprintf('%.2f', deltavbarh) '$'];
annotation('textbox', [xshift + .082 .3 .1 .1], 'String', strparamh, 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 15);

set(gca, 'LineWidth', 3, 'FontSize', 25, 'FontWeight', 'bold');

xlim(xlimvec);
ylim(ylimvec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
choice = 1;
col = [307 498 625 943];
%col = [20 50 80 110];

filename{1} = 'R_1_06.txt';
filename{2} = 'R_1_11.txt';
filename{3} = 'R_1_24.txt';
filename{4} = 'R_1_59_neg.txt';
filename{5} = 'R_1_59_pos.txt';
filename{6} = 'R_2_02_neg.txt';
filename{7} = 'R_2_02_pos.txt';

%Load
data = load(filename{choice});

clr = jet(length(col));
skip = 2;

%subplot(1, 3, 2);
subplot('Position', [xshift + 0.37 .2 0.29 .75]);
hold on;
box on;

nmin = 0;
%[px, py] = plotshaded([-nmin nmin], [3e-4 3e-4; 0.065-3e-4 0.065-3e-4], [0.8 0.8 0.8], '');    

%Get n
switch choice
    case 1
        n = linspace(-240.97, 240.97, 5000).';
        TT = linspace(1.8, 156, 981);
    case 2
        n = linspace(-240.97, 240.97, 5000).';
        TT = linspace(1.8, 156, 981);
    case 3
        n = linspace(-486.3, 488.8, 386).';
        TT = linspace(0, 302, 400);
    case 4
        n = linspace(-658.8, 31.19, 200).';
        TT = linspace(0, 300, 200);
    case 5
        n = linspace(0, 627.6, 182).';
        TT = linspace(0, 300, 200);
    case 6
        n = linspace(-671, 53.5, 150).';
        TT = linspace(0, 300, 200);
    case 7
        n = linspace(-27.1, 671, 130).';
        TT = linspace(0, 300, 110);
end  

switch choice
    case 1
        deg = 1.06;
    case 2
        deg = 1.11;
    case 3
        deg = 1.24;
    case {4, 5}
        deg = 1.59;
    case {6, 7}
        deg = 2.02;
    otherwise
        error('Wrong choice');
end   

for j = 1:length(col)
    T = TT(col(j));
    disp(['T=' num2str(T)]);

    %Extract
    rho = data(:, col(j)) / 25812.80759; %h/e^2
      
    n1 = n(1:skip:end);
    rho1 = rho(1:skip:end);
    index = abs(n1) >= nmin;
    p(j) = plot(n1(index), rho1(index), 's', 'Color', clr(j, :), 'MarkerSize', 9, 'LineWidth', 2, 'DisplayName', strDisp(T)); %'MarkerFaceColor', clr(j, :)
end

plot([0 0], 0.88 * ylimvec, 'k--', 'LineWidth', 3);

xlabel('$n~(10^{10}\mathrm{cm}^{-2})$', 'FontSize', 30, 'Interpreter', 'latex');
%ylabel('$\rho~(h/e^2)$', 'FontSize', 30, 'Interpreter', 'latex');

xticks(xtickvec);
xticklabels(num2cell(xtickvec));

yticks(ytickvec);
yticklabels({});

l = legend(p);
legend boxoff;
set(l, 'Interpreter', 'latex', 'Position', [xshift + .38 .38 .1 .1], 'FontSize', 20);

annotation('textbox', [xshift + .375 .73 .1 .1], 'String', ['$\theta=' num2str(deg) '^\circ$'], 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 22);
annotation('textbox', [xshift + .375 .18 .1 .1], 'String', 'Holes', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 20);
annotation('textbox', [xshift + .582 .18 .1 .1], 'String', 'Electron', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 20);
annotation('textbox', [xshift + .43 .83 .1 .1], 'String', 'Experiment (Ref.~13)', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 22);

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
subplot('Position', [xshift + 0.66 .2 0.29 .75]);
hold on;
box on;

for j = 1:length(Ttheory)
    T = Ttheory(j);
   
	if loadfile
        load(['rho_Planckian_Dirac_deg' num2str(deg) '_T' num2str(T) '_e.mat'], 'T', 'ntheory', 'deg', 'teV', 'Ce', 'vbare', 'rhoDirac_e');
        load(['rho_Planckian_Dirac_deg' num2str(deg) '_T' num2str(T) '_h.mat'], 'T', 'ntheory', 'deg', 'teV', 'Ch', 'vbarh', 'rhoDirac_h');
        
        load(['rho_Planckian_CN_deg' num2str(deg) '_T' num2str(T) '_e.mat'], 'T', 'ntheory', 'deg', 'teV', 'Ce', 'vbare', 'rhoe');
        load(['rho_Planckian_CN_deg' num2str(deg) '_T' num2str(T) '_h.mat'], 'T', 'ntheory', 'deg', 'teV', 'Ch', 'vbarh', 'rhoh');
    else
        %Dirac
        rhoDirac_e = 1 ./ sigmaplanckfun(T, ntheory, deg, teV, Ce, vbare);
        rhoDirac_h = 1 ./ sigmaplanckfun(T, ntheory, deg, teV, Ch, vbarh);
        
        %CN
        rhoe = 1 ./ sigmaPlanckianCN_avg(T, ntheory, deg, teV, Ce, vbare);
        rhoh = 1 ./ sigmaPlanckianCN_avg(T, ntheory, deg, teV, Ch, vbarh);
	end

    if savefile
        save(['rho_Planckian_Dirac_deg' num2str(deg) '_T' num2str(T) '_e.mat'], 'T', 'ntheory', 'deg', 'teV', 'Ce', 'vbare', 'rhoDirac_e');
        save(['rho_Planckian_Dirac_deg' num2str(deg) '_T' num2str(T) '_h.mat'], 'T', 'ntheory', 'deg', 'teV', 'Ch', 'vbarh', 'rhoDirac_h');        
        
        save(['rho_Planckian_CN_deg' num2str(deg) '_T' num2str(T) '_e.mat'], 'T', 'ntheory', 'deg', 'teV', 'Ce', 'vbare', 'rhoe');
        save(['rho_Planckian_CN_deg' num2str(deg) '_T' num2str(T) '_h.mat'], 'T', 'ntheory', 'deg', 'teV', 'Ch', 'vbarh', 'rhoh');
    end    
    
    p(j) = plot([fliplr(-ntheory) ntheory], [fliplr(rhoh) rhoe], '-', 'Color', clr(j, :), 'LineWidth', 3, 'DisplayName', strDisp(T)); %, , 'DisplayName', strDisp{j}
    
    plot([fliplr(-ntheory) ntheory], [fliplr(rhoDirac_h) rhoDirac_e], '--', 'Color', clr(j, :), 'LineWidth', 3);  
end

plot([0 0], 0.88 * ylimvec, 'k--', 'LineWidth', 3); %'Color', [.4 .4 .4]

xticks(xtickvec);
xticklabels(num2cell(xtickvec));

yticks(ytickvec);
yticklabels({});

xlabel('$n~(10^{10}\mathrm{cm}^{-2})$', 'FontSize', 30, 'Interpreter', 'latex');
%ylabel('$\rho~(h/e^2)$', 'FontSize', 30, 'Interpreter', 'latex');

annotation('textbox', [xshift + .665 .18 .1 .1], 'String', 'Holes', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 20);
annotation('textbox', [xshift + .872 .18 .1 .1], 'String', 'Electron', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 20);
annotation('textbox', [0.7487 0.8545 0.1515 0.0755], 'String', 'Planckian theory','Interpreter', 'latex', 'FontSize', 22, 'LineWidth', 2);

strparame{1} = ['$C_e=' sprintf('%.1f', Ce) ' \pm ' sprintf('%.1f', deltaCe) '$'];
strparame{2} = ['$v_F^{(e)}/v_0=' sprintf('%.2f', vbare) ' \pm ' sprintf('%.2f', deltavbare) '$'];
annotation('textbox', [xshift + .825 .28 .1 .1], 'String', strparame, 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 15);
strparamh{1} = ['$C_h=' sprintf('%.2f', Ch) ' \pm ' sprintf('%.2f', deltaCh) '$'];
strparamh{2} = ['$v_F^{(h)}/v_0=' sprintf('%.2f', vbarh) ' \pm ' sprintf('%.2f', deltavbarh) '$'];
annotation('textbox', [xshift + .66 .28 .1 .1], 'String', strparamh, 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 15);

x0 = .09;
dx = 0.29;
y = .81;
annotation('textbox', [0 y .1 .1], 'String', '(a)', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 25);
%annotation('textbox', [x0 + dx y .1 .1], 'String', '(b)', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 25);
%annotation('textbox', [x0 + 2 * dx y .1 .1], 'String', '(c)', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 25);

set(gca, 'LineWidth', 3, 'FontSize', 25, 'FontWeight', 'bold');

xlim(xlimvec);
ylim(ylimvec);

if savefile
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
end

% end
%
% %Shading
%     function [px, py] = plotshaded(x, y, fstr, dispname)
% 
% 	px = [x,fliplr(x)]; % make closed patch
% 	py = [y(1, :), fliplr(y(2, :))];
% 	p(1) = patch(px, py, 1,'FaceColor', fstr, 'EdgeColor', 'none', 'DisplayName', dispname); %'FaceAlpha', 0.2
%     end