%extract data vs T

clear;

savepdf = 1;
loadfile = 1;

% choice = 3;
%row = 193:197; %193:5:235

%row = 193:197;
% row = 193;
% 
% filename{1} = 'R_1_06.txt';
% filename{2} = 'R_1_11.txt';
% filename{3} = 'R_1_24.txt';
% filename{4} = 'R_1_59_neg.txt';
% filename{5} = 'R_1_59_pos.txt';
% filename{6} = 'R_2_02_neg.txt';
% filename{7} = 'R_2_02_pos.txt';
% 
% %Load
% data = load(filename{choice});
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
% switch choice
%     case 1
%         deg = 1.06;
%     case 2
%         deg = 1.11;
%     case 3
%         deg = 1.24;
%     case {4, 5}
%         deg = 1.59;
%     case {6, 7}
%         deg = 2.02;
%     otherwise
%         error('Wrong choice');
% end 
% 
% % figure();
% % hold on;
% % box on;
% 
% %clr = 'kbrmc';
% clr = jet(length(row));

deg = 1.24;

%Theory
Ttheory = [linspace(1, 20, 40) linspace(20.5, 150, 60)]; %linspace(1, 50, 120)

sigmas = 104.35 * 1.0; 
nimp = 5.7; %21.39
nrms = 3.61;
alpha0 = 0.7;
betaAeV = 3.6; %3.6

degstr = num2str(deg);
degstr = strrep(degstr, '.', '_');
% if savefile, print('-dpdf', ['rho_deg' degstr '_pos_vsT_exp.pdf']); end

figure();
hold on;
box on;

% for j = 1:length(row)
%     n = nn(row(j))
j = 1;
n = -0.02;

    if loadfile
        load(['rho_tot_emt_deg' num2str(deg) '_n' sprintf('%.2f', n) '_CN.mat'], 'Ttheory', 'n', 'deg', 'sigmas', 'nimp', 'nrms', 'alpha0', 'betaAeV', 'rhotot');
%        load(['rho_tot_extend_deg' num2str(deg) '_n' sprintf('%.2f', n) '_CN.mat'], 'Ttheory_extend', 'n', 'deg', 'sigmas', 'nimp', 'nrms', 'alpha0', 'betaAeV', 'rhotot_extend');
    else
        %rhotot = 1 ./ sigmatotpolT0sigmas_emt_nrms(Ttheory, deg, n, sigmas, nimp, nrms, alpha0, betaAeV);
        rhotot = 1 ./ sigmatotemtsigmas_nrms(Ttheory, deg, n, sigmas, nimp, nrms, alpha0, betaAeV);
        
        save(['rho_tot_emt_deg' num2str(deg) '_n' sprintf('%.2f', n) '_CN.mat'], 'Ttheory', 'n', 'deg', 'sigmas', 'nimp', 'nrms', 'alpha0', 'betaAeV', 'rhotot')
    end
        
    if sprintf('%.0f', abs(n)) == '0'
        txtunit = '';
    else
        txtunit = '\times 10^{10}\mathrm{cm}^{-2}';
    end    
    
    %clr(j, :)
    p(j) = plot(Ttheory, rhotot, '-', 'Color', [0 0 0], 'LineWidth', 3);
%    p(j) = plot([linspace(1, 20, 40) linspace(20.5, 150, 60)], rhotot, '-', 'Color', [0 0 0], 'LineWidth', 3);
% end

T1 = 3.5;
T2 = 16;
T3 = 65;
ybound = [0.0 0.145];
plot(T1 * [1 1], ybound, 'k--', 'LineWidth', 2);
plot(T2 * [1 1], ybound, 'k--', 'LineWidth', 2);
plot(T3 * [1 1], ybound, 'k--', 'LineWidth', 2);

xlabel('$T~(K)$', 'FontSize', 30, 'Interpreter', 'latex');
ylabel('$\rho~(h/e^2)$', 'FontSize', 30, 'Interpreter', 'latex');

% l = legend(p);
% legend boxoff;
% set(l, 'Interpreter', 'latex', 'Location', 'NorthEast', 'FontSize', 20);       
        
% annotation('textbox', [0.2 0.8133 0.1227 0.0867], 'String', ['$\theta=' num2str(deg) '^\circ$'], 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 22);
% 
% txtparam{1} = ['$n_\mathrm{imp}=' num2str(nimp) '\times 10^{10}\mathrm{cm}^{-2}$'];
% txtparam{2} = ['$\sigma_s=' num2str(sigmas) '~e^2/h$'];
% txtparam{3} = ['$n_\mathrm{rms}=' num2str(nrms) '\times 10^{10}\mathrm{cm}^{-2}$'];
% 
% annotation('textbox', [0.56 0.27 0.1227 0.0867], 'String', txtparam, 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 18);
% 

% xtip=0.4; ytip=0.13;   % arrow tip coordinates (normalized units)
% w=0.1;      %box width
% h=0.1;      %box height
% offset=0.1;
% str={'$T > T_F$','impurity'};
% 
% x = [xtip-offset xtip];     % arrows start and end coordinates
% y = [ytip+offset ytip];     % I've just offset by 0.1 in x and y. 
% a=annotation('textarrow',x,y,'Color','black');
% b=annotation('textbox',[xtip-w-0.1 ytip+0.1 w h],'String',str, 'Color','black', 'EdgeColor','black', 'Interpreter', 'latex','FontSize', 22);

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

if savepdf, print('-dpdf', 'regime.pdf'); end