%function plotFig_1pt11()

clear;

savefile = 1;
loadfile = 1;

deg = 1.11;

% scrsz = get(0, 'ScreenSize');
% figure('Position', [1 1 .6 * scrsz(3) .5 * scrsz(4)]);
figure('Position', [1 1 .6 * 1920 .35 * 1080]); %[1 1 .6 * 1920 .5 * 1080]

% strDisp{1} = '$|n|=10^{11} \mathrm{cm}^{-2}$';
% strDisp{2} = '$|n|=2\times 10^{11} \mathrm{cm}^{-2}$';
% strDisp{3} = '$|n|=5\times 10^{11} \mathrm{cm}^{-2}$';
% strDisp{4} = '$|n|=10^{12} \mathrm{cm}^{-2}$';
strDisp{1} = '$n=-10^{11} \mathrm{cm}^{-2}$';
strDisp{2} = '$n=-2\times 10^{11} \mathrm{cm}^{-2}$';
strDisp{3} = '$n=-3\times 10^{11} \mathrm{cm}^{-2}$';
strDisp{4} = '$n=-4\times 10^{11} \mathrm{cm}^{-2}$';

%xtickvec = [-100 -50 0 50 100];
xtickvec = [0 50 100];
%xtickvecmid = [-200 -100 0 100 200];
ytickvec = 0.1:0.1:0.5;

%xlimvec = [-160 160];
xlimvec = [0 160];
ylimvec = [0 0.5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Theory
%nn = [10 20 50 100];
ntheory = [10 20 30 40];
Ttheory = linspace(.1, 160, 30);
teV = 3.4;
% betaAeV = 6.366;
% vbar = 0.117;

%Theory e-ph
% load(['parameph_deg' num2str(deg) '_e.mat']);
% paramcell = num2cell(param);
% paramdeltacell = num2cell(paramdelta);
% [betaAtildeeVe, vbare] = paramcell{:};
% [deltabetaAtildeeVe, deltavbare] = paramdeltacell{:};

load(['parameph_deg' num2str(deg) '_h.mat']);
paramcell = num2cell(param);
paramdeltacell = num2cell(paramdelta);
[betaAtildeeVh, vbarh] = paramcell{:};
[deltabetaAtildeeVh, deltavbarh] = paramdeltacell{:};

clr = jet(numel(ntheory));

%subplot(1, 3, 1);
subplot('Position', [0.08 .2 0.29 .75]); %[0.08 .13 0.29 .8]
hold on;
box on;

for j = 1:4 %coz Dirac model is weakly density dependent anyway
    n = ntheory(j);

    if loadfile
%         load(['rhoeph_vsT_deg' num2str(deg) '_n' num2str(n) '_Dirac_e.mat'], 'Ttheory', 'n', 'deg', 'betaAtildeeVe', 'vbare', 'rhoeph_Dirac_e');
        load(['rhoeph_vsT_deg' num2str(deg) '_n' num2str(n) '_Dirac_h.mat'], 'Ttheory', 'n', 'deg', 'betaAtildeeVh', 'vbarh', 'rhoeph_Dirac_h');
%         load(['rhoeph_vsT_deg' num2str(deg) '_n' num2str(n) '_CN_e.mat'], 'Ttheory', 'n', 'deg', 'betaAtildeeVe', 'vbare', 'rhoeph_CN_e');
        load(['rhoeph_vsT_deg' num2str(deg) '_n' num2str(n) '_CN_h.mat'], 'Ttheory', 'n', 'deg', 'betaAtildeeVh', 'vbarh', 'rhoeph_CN_h');
    else
%         rhoeph_Dirac_e = 1 ./ sigmaephfun(Ttheory, n, deg, teV, betaAtildeeVe, vbare);
        rhoeph_Dirac_h = 1 ./ sigmaephfun(Ttheory, n, deg, teV, betaAtildeeVh, vbarh);
        %rhoeph_Dirac_BG = rhoephfun(Ttheory, n, deg, betaAeV);
        %rhotot_Dirac = 1 ./ sigmatotemt_nrms(Ttheory, deg, n, sigmas, nimp, nrms, alpha0, betaAeV);
        
%         rhoeph_CN_e = 1./ sigmaephCN_exact_avg_iso(Ttheory, n, deg, [], teV, betaAtildeeVe, vbare);
        rhoeph_CN_h = 1./ sigmaephCN_exact_avg_iso(Ttheory, n, deg, [], teV, betaAtildeeVh, vbarh);
    end
    
    if savefile
%         save(['rhoeph_vsT_deg' num2str(deg) '_n' num2str(n) '_Dirac_e.mat'], 'Ttheory', 'n', 'deg', 'betaAtildeeVe', 'vbare', 'rhoeph_Dirac_e');
        save(['rhoeph_vsT_deg' num2str(deg) '_n' num2str(n) '_Dirac_h.mat'], 'Ttheory', 'n', 'deg', 'betaAtildeeVh', 'vbarh', 'rhoeph_Dirac_h');
%         save(['rhoeph_vsT_deg' num2str(deg) '_n' num2str(n) '_CN_e.mat'], 'Ttheory', 'n', 'deg', 'betaAtildeeVe', 'vbare', 'rhoeph_CN_e');
        save(['rhoeph_vsT_deg' num2str(deg) '_n' num2str(n) '_CN_h.mat'], 'Ttheory', 'n', 'deg', 'betaAtildeeVh', 'vbarh', 'rhoeph_CN_h');
    end
    
    %p(j) = plot([fliplr(-Ttheory) NaN Ttheory], [fliplr(rhoeph_Dirac_h) NaN rhoeph_Dirac_e], '--', 'Color', clr(j, :), 'LineWidth', 3, 'DisplayName', sprintf('Electron-phonon\nDirac model')); %['$n=' sprintf('%.0f', abs(n)) '\times 10^{10}\mathrm{cm}^{-2}$']
    index = rhoeph_Dirac_h < 0.417;
    p(j) = plot(Ttheory(index), rhoeph_Dirac_h(index), '--', 'Color', clr(j, :), 'LineWidth', 3, 'DisplayName', sprintf('Electron-phonon\nDirac model')); %['$n=' sprintf('%.0f', abs(n)) '\times 10^{10}\mathrm{cm}^{-2}$']
    
%     rhoeph_Dirac_BG = rhoephfun(Ttheory, n, deg, betaAeV);
%     plot([fliplr(-Ttheory) NaN Ttheory], [fliplr(rhoeph_Dirac_BG) NaN rhoeph_Dirac_BG], '-.', 'Color', clr(j, :), 'LineWidth', 3);
    
    %plot([fliplr(-Ttheory) NaN Ttheory], [fliplr(rhoeph_CN_h) NaN rhoeph_CN_e], '-', 'Color', clr(j, :), 'LineWidth', 3, 'DisplayName', sprintf('Electron-phonon\nCastro-Neto model'));
    plot(Ttheory, rhoeph_CN_h, '-', 'Color', clr(j, :), 'LineWidth', 3, 'DisplayName', sprintf('Electron-phonon\nCastro-Neto model'));
end

%plot([0 0], ylimvec, 'k--', 'LineWidth', 3);

xticks(xtickvec);
%xticklabels({'200','100','0','100','200'});
xticklabels(num2cell(abs(xtickvec)));

xlabel('$T$~(K)', 'FontSize', 30, 'Interpreter', 'latex');
ylabel('$\rho~(h/e^2)$', 'FontSize', 30, 'Interpreter', 'latex');

% annotation('textbox', [.085 .75 .1 .1], 'String', ['$\theta=' num2str(deg) '^\circ$'], 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 22);
% annotation('textbox', [.085 .1 .1 .1], 'String', 'Holes', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 20);
% annotation('textbox', [.292 .1 .1 .1], 'String', 'Electron', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 20);
annotation('textbox', [0.1578 0.8439 0.1329 0.0761], 'String', 'Phonon theory', 'Interpreter', 'latex', 'FontSize', 22, 'LineWidth', 2); %, 'LineStyle', 'none'
% annotation('textbox', [.085 .1 .1 .1], 'String', {'Dashed lines:','Dirac model'}, 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 15);

% rad = (deg / 180) * pi;
% gamma = 1 ./ (2 * tan(rad/2));
% betaAtildeeVe = betaAeVe * vbare * gamma;
% deltabetaAtildeeVe = (betaAeVe * deltavbare + deltabetaAeVe * vbare) * gamma;
% betaAtildeeVh = betaAeVh * vbarh * gamma;
% deltabetaAtildeeVh = (betaAeVh * deltavbarh + deltabetaAeVh * vbarh) * gamma;

% strparame{1} = ['$\tilde{\beta}_A^{(e)}=' sprintf('%.1f', betaAtildeeVe) ' \pm ' sprintf('%.1f', deltabetaAtildeeVe) '$~eV'];
% strparame{2} = ['$v_F^{(e)}/v_0=' sprintf('%.3f', vbare) ' \pm ' sprintf('%.3f', deltavbare) '$'];
% annotation('textbox', [.23 .75 .1 .1], 'String', strparame, 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 15);
strparamh{1} = ['$\tilde{\beta}_A^{(h)}=' sprintf('%.0f', betaAtildeeVh) ' \pm ' sprintf('%.0f', deltabetaAtildeeVh) '$~eV'];
strparamh{2} = ['$v_F^{(h)}/v_0=' sprintf('%.2f', vbarh) ' \pm ' sprintf('%.2f', deltavbarh) '$'];
annotation('textbox', [0.22 0.25  .1 .1], 'String', strparamh, 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 17);

set(gca, 'LineWidth', 3, 'FontSize', 25, 'FontWeight', 'bold');

xlim(xlimvec);
ylim(ylimvec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
choice = 2;
% row{1} = [2604 2708 3019 3538];
% row{2} = fliplr([1463 1982 2293 2397]);
row{1} = [2604 2708 2812 2915];
row{2} = fliplr([2086 2189 2293 2397]);

filename{1} = 'R_1_06.txt';
filename{2} = 'R_1_11.txt';
filename{3} = 'R_1_24.txt';
filename{4} = 'R_1_59_neg.txt';
filename{5} = 'R_1_59_pos.txt';
filename{6} = 'R_2_02_neg.txt';
filename{7} = 'R_2_02_pos.txt';

%Load
data = load(filename{choice});

clr = jet(length(row{1}));
skipvec = [2 2];

%subplot(1, 3, 2);
subplot('Position', [0.37 .2 0.29 .75]); %[0.37 .13 0.29 .8]
hold on;
box on;

Tmin = 0;
%[px, py] = plotshaded([-Tmin Tmin], [3e-4 3e-4; 0.065-3e-4 0.065-3e-4], [0.8 0.8 0.8], '');    

%Get n
switch choice
    case 1
        nn = linspace(-240.97, 240.97, 5000).';
        T = linspace(1.8, 156, 981);
    case 2
        nn = linspace(-240.97, 240.97, 5000).';
        T = linspace(1.8, 156, 981);
    case 3
        nn = linspace(-486.3, 488.8, 386).';
        T = linspace(0, 302, 400);
    case 4
        nn = linspace(-658.8, 31.19, 200).';
        T = linspace(0, 300, 200);
    case 5
        nn = linspace(0, 627.6, 182).';
        T = linspace(0, 300, 200);
    case 6
        nn = linspace(-671, 53.5, 150).';
        T = linspace(0, 300, 200);
    case 7
        nn = linspace(-27.1, 671, 130).';
        T = linspace(0, 300, 110);
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
    
for k = 2:2 %1:2
    cur_row = row{k};
    
    for j = 1:length(cur_row)
        n = nn(cur_row(j));
        disp(['n=' num2str(n)]);

        %Extract
        rho = data(cur_row(j), :) / 25812.80759; %h/e^2

        skip = skipvec(k);
        %p(j) = plot(T, rho, clr(j), 'LineWidth', 2, 'DisplayName', ['$n=' sprintf('%.1f', n) '\times 10^{10}\mathrm{cm}^{-2}$']);
        %strDisp = ['$|n|=' sprintf('%.0f', n) '\times 10^{10}\mathrm{cm}^{-2}$'];

        T1 = T(1:skip:end);
        rho1 = rho(1:skip:end);
        index = T1 >= Tmin;
        %p(j) = plot(T1(index) * (-1) ^ (k + 1), rho1(index), 's', 'Color', clr(j, :), 'MarkerSize', 9, 'LineWidth', 2, 'DisplayName', strDisp{j}); %'MarkerFaceColor', clr(j, :)
        p(j) = plot(T1(index), rho1(index), 's', 'Color', clr(j, :), 'MarkerSize', 9, 'LineWidth', 2, 'DisplayName', strDisp{j}); %'MarkerFaceColor', clr(j, :)
    end
end

%Theory
%Ttheory_eph_Dirac = linspace(Tmin, 300, 40);
% Ttheory_Planckian_Dirac = linspace(Tmin, 300, 40);
% betaAeV = 3.6;
% 
% %Get n
% switch choice
%     case 1
%         nn = linspace(-240.97, 240.97, 5000).';
%         T = linspace(1.8, 156, 981);
%     case 2
%         nn = linspace(-240.97, 240.97, 5000).';
%         T = linspace(1.8, 156, 981);
%     case 3
%         nn = linspace(-486.3, 488.8, 386).';
%         T = linspace(0, 302, 400);
%     case 4
%         nn = linspace(-658.8, 31.19, 200).';
%         T = linspace(0, 300, 200);
%     case 5
%         nn = linspace(0, 627.6, 182).';
%         T = linspace(0, 300, 200);
%     case 6
%         nn = linspace(-671, 53.5, 150).';
%         T = linspace(0, 300, 200);
%     case 7
%         nn = linspace(-27.1, 671, 130).';
%         T = linspace(0, 300, 110);
% end
%
% nn = [10 20 50 100];
% 
% for j = 1:4 %coz Dirac model is weakly density dependent anyway
%     cur_row = row{k};
%     %n = nn(cur_row(j));
%     n = nn(j);
% 
%     if loadfile
%         load(['rhoeph_vsT_deg' num2str(deg) '_n' num2str(n) '_Dirac.mat'], 'Ttheory_eph_Dirac', 'n', 'deg', 'betaAeV', 'vbar', 'rhoeph_Dirac');    
%     end    
%     %rhoeph = 1 ./ arrayfun(@(TT) sigmaephfun(TT, n, deg, betaAeV, vbar), Ttheory_eph_Dirac);
%     %rhoeph = rhoephfun(Ttheory_eph_Dirac, n, deg, betaAeV);
%     %rhotot = 1 ./ sigmatotemt_nrms(Ttheory, deg, n, sigmas, nimp, nrms, alpha0, betaAeV);
%     index = Ttheory_eph_Dirac >= Tmin;
%     p(5) = plot([fliplr(-Ttheory_eph_Dirac(index)) NaN Ttheory_eph_Dirac(index)], [fliplr(rhoeph_Dirac(index)) NaN rhoeph_Dirac(index)], '-', 'Color', 'k', 'LineWidth', 3, 'DisplayName', sprintf('Electron-phonon\nDirac model')); %['$n=' sprintf('%.0f', abs(n)) '\times 10^{10}\mathrm{cm}^{-2}$']
%     
% %     rhoPlanckianDirac = 1 ./ sigmaplanckfun(Ttheory_Planckian_Dirac, n, deg, vbar);
% %     plot([fliplr(-Ttheory_Planckian_Dirac) NaN Ttheory_Planckian_Dirac], [fliplr(rhoPlanckianDirac) NaN rhoPlanckianDirac], '--', 'Color', clr(j, :), 'LineWidth', 3);
% end

%plot([0 0], ylimvec, 'k--', 'LineWidth', 3);

xlabel('$T$~(K)', 'FontSize', 30, 'Interpreter', 'latex');
%ylabel('$\rho~(h/e^2)$', 'FontSize', 30, 'Interpreter', 'latex');

xticks(xtickvec);
%xticklabels({'200','100','0','100','200'});
xticklabels(num2cell(abs(xtickvec)));

yticks(ytickvec);
yticklabels({});

l = legend(p);
%l = legend('$10^{11} \mathrm{cm}^{-2}$', '$2\times 10^{11} \mathrm{cm}^{-2}$', '$5\times 10^{11} \mathrm{cm}^{-2}$', '$10^{12} \mathrm{cm}^{-2}$', '', 'Dirac theory');
legend boxoff;
set(l, 'Interpreter', 'latex', 'Location', 'SouthEast', 'FontSize', 17);
%set(l, 'Interpreter', 'latex', 'Position', [.385 .25 .1 .1], 'FontSize', 15);

annotation('textbox', [.38 .73 .1 .1], 'String', ['$\theta=' num2str(deg) '^\circ$'], 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 22);
% annotation('textbox', [.375 .1 .1 .1], 'String', 'Holes', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 20);
% annotation('textbox', [.582 .1 .1 .1], 'String', 'Electron', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 20);
annotation('textbox', [0.4248 0.8307 0.1889 0.0868], 'String', 'Experiment (Ref.~13)', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 22);
%annotation('textbox', [.53 .26 .1 .1], 'String', {'Dashed lines:', 'Planckian theory','Dirac model'}, 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 15);

set(gca, 'LineWidth', 3, 'FontSize', 25, 'FontWeight', 'bold');

xlim(xlimvec);
ylim(ylimvec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Theory Planckian
%vbar = 0.117;
% load(['paramPlanckian_deg' num2str(deg) '_e.mat']);
% paramcell = num2cell(param);
% paramdeltacell = num2cell(paramdelta);
% [Ce, vbare] = paramcell{:};
% [deltaCe, deltavbare] = paramdeltacell{:};

load(['paramPlanckian_deg' num2str(deg) '_h.mat']);
paramcell = num2cell(param);
paramdeltacell = num2cell(paramdelta);
[Ch, vbarh] = paramcell{:};
[deltaCh, deltavbarh] = paramdeltacell{:};

clr = jet(length(ntheory));

%subplot(1, 3, 3);
subplot('Position', [0.66 .2 0.29 .75]); %[0.66 .13 0.29 .8]
hold on;
box on;

% p = [];
% yyaxis left;
% ytickvec = linspace(0, .2, 5);
% yticks(ytickvec);
% yticklabels({});

%ylim([0 .6]);

% yyaxis right;

%plot([fliplr(-Ttheory_eph_Dirac) NaN Ttheory_eph_Dirac], [fliplr(rhoeph_Dirac) NaN rhoeph_Dirac], '--', 'Color', clr(1, :), 'LineWidth', 3);

for j = 1:length(ntheory)
    n = ntheory(j);
   
	if loadfile
%         load(['rho_Planckian_Dirac_deg' num2str(deg) '_n' sprintf('%.2f', n) '_e.mat'], 'Ttheory', 'n', 'deg', 'teV', 'Ce', 'vbare', 'rhoDirac_e');
        load(['rho_Planckian_Dirac_deg' num2str(deg) '_n' sprintf('%.2f', n) '_h.mat'], 'Ttheory', 'n', 'deg', 'teV', 'Ch', 'vbarh', 'rhoDirac_h');
        
%         load(['rho_Planckian_CN_deg' num2str(deg) '_n' sprintf('%.2f', n) '_e.mat'], 'Ttheory', 'n', 'deg', 'teV', 'Ce', 'vbare', 'rhoe');
        load(['rho_Planckian_CN_deg' num2str(deg) '_n' sprintf('%.2f', n) '_h.mat'], 'Ttheory', 'n', 'deg', 'teV', 'Ch', 'vbarh', 'rhoh');
    else
        %Dirac
%         rhoDirac_e = 1 ./ sigmaplanckfun(Ttheory, n, deg, teV, Ce, vbare);
        rhoDirac_h = 1 ./ sigmaplanckfun(Ttheory, n, deg, teV, Ch, vbarh);
        
        %CN
%         rhoe = 1 ./ sigmaPlanckianCN_avg(Ttheory, n, deg, teV, Ce, vbare);
        rhoh = 1 ./ sigmaPlanckianCN_avg(Ttheory, n, deg, teV, Ch, vbarh);
	end

    if savefile
%         save(['rho_Planckian_Dirac_deg' num2str(deg) '_n' sprintf('%.2f', n) '_e.mat'], 'Ttheory', 'n', 'deg', 'teV', 'Ce', 'vbare', 'rhoDirac_e');
        save(['rho_Planckian_Dirac_deg' num2str(deg) '_n' sprintf('%.2f', n) '_h.mat'], 'Ttheory', 'n', 'deg', 'teV', 'Ch', 'vbarh', 'rhoDirac_h');        
        
%         save(['rho_Planckian_CN_deg' num2str(deg) '_n' sprintf('%.2f', n) '_e.mat'], 'Ttheory', 'n', 'deg', 'teV', 'Ce', 'vbare', 'rhoe');
        save(['rho_Planckian_CN_deg' num2str(deg) '_n' sprintf('%.2f', n) '_h.mat'], 'Ttheory', 'n', 'deg', 'teV', 'Ch', 'vbarh', 'rhoh');
    end    
    
    %p(j) = plot([fliplr(-Ttheory) Ttheory], [fliplr(rhoh) rhoe], '-', 'Color', clr(j, :), 'LineWidth', 3, 'DisplayName', strDisp{j}); %, , 'DisplayName', strDisp{j}
    p(j) = plot(Ttheory, rhoh, '-', 'Color', clr(j, :), 'LineWidth', 3, 'DisplayName', strDisp{j}); %, , 'DisplayName', strDisp{j}
    
    %plot([fliplr(-Ttheory) Ttheory], [fliplr(rhoDirac_h) rhoDirac_e], '--', 'Color', clr(j, :), 'LineWidth', 3);
    plot(Ttheory, rhoDirac_h, '--', 'Color', clr(j, :), 'LineWidth', 3);
end

%plot([0 0], ylimvec, 'k--', 'LineWidth', 3); %'Color', [.4 .4 .4]

% xticks([-40 -20 0 20 40]);
% %xticklabels({'200','100','0','100','200'});
% xticklabels(num2cell(abs([-40 -20 0 20 40])));

xticks(xtickvec);
xticklabels(num2cell(abs(xtickvec)));

yticks(ytickvec);
%yticklabels(num2cell(abs(ytickvec)));
yticklabels({});

% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.YAxis(2).Color = 'k';

xlabel('$T$~(K)', 'FontSize', 30, 'Interpreter', 'latex');
%ylabel('$\rho~(h/e^2)$', 'FontSize', 30, 'Interpreter', 'latex');

% l = legend(p);
% %l = legend('$10^{11} \mathrm{cm}^{-2}$', '$2\times 10^{11} \mathrm{cm}^{-2}$', '$5\times 10^{11} \mathrm{cm}^{-2}$', '$10^{12} \mathrm{cm}^{-2}$', '', 'Dirac theory');
% legend boxoff;
% set(l, 'Interpreter', 'latex', 'Position', [.69 .225 .1 .1], 'FontSize', 15); %'Location', 'SouthWest'

% annotation('textbox', [.67 .75 .1 .1], 'String', ['$\theta=' num2str(deg) '^\circ$'], 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 22);
% annotation('textbox', [.665 .1 .1 .1], 'String', 'Holes', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 20);
% annotation('textbox', [.872 .1 .1 .1], 'String', 'Electron', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 20);
annotation('textbox', [0.7313 0.8439 0.1515 0.0761], 'String', 'Planckian theory', 'Interpreter', 'latex', 'FontSize', 22, 'LineWidth', 2); %, 'LineStyle', 'none'
% annotation('textbox', [.665 .11 .1 .1], 'String', {'Dashed lines:','Dirac model'}, 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 15);

% strparame{1} = ['$C_e=' sprintf('%.1f', Ce) ' \pm ' sprintf('%.1f', deltaCe) '$'];
% strparame{2} = ['$v_F^{(e)}/v_0=' sprintf('%.3f', vbare) ' \pm ' sprintf('%.3f', deltavbare) '$'];
% annotation('textbox', [.81 .75 .1 .1], 'String', strparame, 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 15);
strparamh{1} = ['$C_h=' sprintf('%.1f', Ch) ' \pm ' sprintf('%.1f', deltaCh) '$'];
strparamh{2} = ['$v_F^{(h)}/v_0=' sprintf('%.2f', vbarh) ' \pm ' sprintf('%.2f', deltavbarh) '$'];
annotation('textbox', [0.8 0.22 .1 .1], 'String', strparamh, 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 17);

x0 = .09;
dx = 0.29;
y = .81;
annotation('textbox', [0 y .1 .1], 'String', '(b)', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 25);
%annotation('textbox', [x0 + dx y .1 .1], 'String', '(b)', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 25);
%annotation('textbox', [x0 + 2 * dx y .1 .1], 'String', '(c)', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 25);

set(gca, 'LineWidth', 3, 'FontSize', 25, 'FontWeight', 'bold');

%xlim([-160 160]);
xlim(xlimvec);
ylim(ylimvec);

% degstr = num2str(deg);
% degstr = strrep(degstr, '.', '_');
%if savefile, print('-dpdf', ['rho_vsT_beyondDirac_deg' num2str(deg) '.pdf']); end
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
    print('-dpdf', ['rho_vsT_beyondDirac_deg' num2str(deg) '.pdf']);

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