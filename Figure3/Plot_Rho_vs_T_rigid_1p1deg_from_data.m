clc
clear all 

load Rhorigid_vs_T_1p1deg_new1.mat 

figure;
box on;
hold on;

xlabel('$T(K)$', 'FontSize', 30, 'Interpreter', 'latex');
ylabel('$\rho_{e-ph}(h/e^2)$', 'FontSize', 30, 'Interpreter', 'latex');

pl1 = plot(TT, rhoephns(1,:), 'r-', 'LineWidth', 4);
pl2 = plot(TT, rhoephns(2,:), 'b-', 'LineWidth', 4);
pl3 = plot(TT, rhoephns(3,:),'-', 'Color','[.2 .7 .2]', 'LineWidth', 4);

plot([Tf(1)/2 Tf(1)/2],[0 1],'--r', 'LineWidth', 2)
plot([Tf(2)/2 Tf(2)/2],[0 1],'--b', 'LineWidth', 2)
plot([Tf(3)/2 Tf(3)/2],[0 1],'--','Color','[.2 .7 .2]', 'LineWidth', 2)

plot([Tbg(1)/4 Tbg(1)/4],[0 1],':r', 'LineWidth', 2)
plot([Tbg(2)/4 Tbg(2)/4],[0 1],':b', 'LineWidth', 2)
plot([Tbg(3)/4 Tbg(3)/4],[0 1],':','Color','[.2 .7 .2]', 'LineWidth', 2)

h = gca;
h.XMinorTick='on';
h.YMinorTick='on';

limm=1.01;
xlim([0 max(TT)/limm]);
ylim([-0 max(max(rhoephns))/(limm)]);

legend([pl1 pl2 pl3],'$10^{10}$cm$^{-2}$','$5\times 10^{10}$cm$^{-2}$','$10\times 10^{10}$cm$^{-2}$', 'Interpreter', 'latex','Location','SouthEast');
legend boxoff

set(gca, 'LineWidth', 3, 'FontSize', 30, 'FontWeight', 'bold');

strx = 0.1;
stry = 0.7;
cc = 30;

strcolor = 'Red';
dim = [strx+0.025 stry-0.15 0.3 0.3];
str = '$\frac{T_{F}}{2}$';
annotation('textbox',dim,'String',str,'LineStyle','none','FontSize',cc-5,'FontWeight', 'bold','Interpreter', 'latex','FontName','Times New Roman', 'Color', strcolor);

strcolor = 'Blue';
dim = [strx+0.1 stry-0.15 0.3 0.3];
str = '$\frac{T_{F}}{2}$';
annotation('textbox',dim,'String',str,'LineStyle','none','FontSize',cc-5,'FontWeight', 'bold','Interpreter', 'latex','FontName','Times New Roman', 'Color', strcolor);

strcolor = '[.2 .7 .2]';
dim = [strx+0.17 stry-0.15 0.3 0.3];
str = '$\frac{T_{F}}{2}$';
annotation('textbox',dim,'String',str,'LineStyle','none','FontSize',cc-5,'FontWeight', 'bold','Interpreter', 'latex','FontName','Times New Roman', 'Color', strcolor);


strcolor = 'Red';
dim = [strx+0.25 stry-0.45 0.3 0.3];
str = '$\frac{T_{BG}}{4}$';
annotation('textbox',dim,'String',str,'LineStyle','none','FontSize',cc-5,'FontWeight', 'bold','Interpreter', 'latex','FontName','Times New Roman', 'Color', strcolor);

strcolor = 'Blue';
dim = [strx+0.48 stry-0.45 0.3 0.3];
str = '$\frac{T_{BG}}{4}$';
annotation('textbox',dim,'String',str,'LineStyle','none','FontSize',cc-5,'FontWeight', 'bold','Interpreter', 'latex','FontName','Times New Roman', 'Color', strcolor);

strcolor = '[.2 .7 .2]';
dim = [strx+0.7 stry-0.45 0.3 0.3];
str = '$\frac{T_{BG}}{4}$';
annotation('textbox',dim,'String',str,'LineStyle','none','FontSize',cc-5,'FontWeight', 'bold','Interpreter', 'latex','FontName','Times New Roman', 'Color', strcolor);


strcolor = 'Black';
dim = [strx+0.3 stry-0.25 0.3 0.3];
str = '$\theta\sim 1.1^\circ$';
annotation('textbox',dim,'String',str,'LineStyle','none','FontSize',cc-5,'FontWeight', 'bold','Interpreter', 'latex','FontName','Times New Roman', 'Color', strcolor);

strcolor = 'Black';
dim = [strx+0.3 stry-0.32 0.3 0.3];
str = '\textit{Interband}';
annotation('textbox',dim,'String',str,'LineStyle','none','FontSize',cc-5,'FontWeight', 'bold','Interpreter', 'latex','FontName','Times New Roman', 'Color', strcolor);

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
print('-dpdf', 'rho_vs_T_interband.pdf');

% Restore the previous settings
set(gcf,'PaperType',prePaperType);
set(gcf,'PaperUnits',prePaperUnits);
set(gcf,'Units',preUnits);
set(gcf,'PaperPosition',prePaperPosition);
set(gcf,'PaperSize',prePaperSize); 
