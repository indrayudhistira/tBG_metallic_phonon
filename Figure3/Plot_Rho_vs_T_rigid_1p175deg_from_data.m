%%
clc
clear all 

load Rhorigid_vs_T_1p75deg_new1.mat 

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

limm=11.5;
xlim([0 max(TT)/limm]);
ylim([-0 max(max(rhoephns))/(limm)]);

legend([pl1 pl2 pl3],'$10^{12}$cm$^{-2}$','$5\times 10^{12}$cm$^{-2}$','$10\times 10^{12}$cm$^{-2}$', 'Interpreter', 'latex','Location','SouthEast');
legend boxoff

set(gca, 'LineWidth', 3, 'FontSize', 30, 'FontWeight', 'bold');

strx = 0.19;
stry = 0.7;
cc = 30;

strcolor = 'Red';
dim = [strx+0.025 stry-0.15 0.3 0.3];
str = '$\frac{T_{BG}}{4}$';
annotation('textbox',dim,'String',str,'LineStyle','none','FontSize',cc-5,'FontWeight', 'bold','Interpreter', 'latex','FontName','Times New Roman', 'Color', strcolor);

strcolor = 'Blue';
dim = [strx+0.09 stry-0.15 0.3 0.3];
str = '$\frac{T_{BG}}{4}$';
annotation('textbox',dim,'String',str,'LineStyle','none','FontSize',cc-5,'FontWeight', 'bold','Interpreter', 'latex','FontName','Times New Roman', 'Color', strcolor);

strcolor = '[.2 .7 .2]';
dim = [strx+0.15 stry-0.15 0.3 0.3];
str = '$\frac{T_{BG}}{4}$';
annotation('textbox',dim,'String',str,'LineStyle','none','FontSize',cc-5,'FontWeight', 'bold','Interpreter', 'latex','FontName','Times New Roman', 'Color', strcolor);


strcolor = 'Red';
dim = [strx+0.6 stry-0.45 0.3 0.3];
str = '$\frac{T_{F}}{2}$';
annotation('textbox',dim,'String',str,'LineStyle','none','FontSize',cc-5,'FontWeight', 'bold','Interpreter', 'latex','FontName','Times New Roman', 'Color', strcolor);

strcolor = 'Blue';
dim = [strx+0.27 stry-0.1 0.3 0.3];
str = '$T_F\sim 45 K$';
annotation('textbox',dim,'String',str,'LineStyle','none','FontSize',cc-10,'FontWeight', 'bold','Interpreter', 'latex','FontName','Times New Roman', 'Color', strcolor);

strcolor = '[.2 .7 .2]';
dim = [strx+0.47 stry-0.1 0.3 0.3];
str = '$T_F\sim 64K$';
annotation('textbox',dim,'String',str,'LineStyle','none','FontSize',cc-10,'FontWeight', 'bold','Interpreter', 'latex','FontName','Times New Roman', 'Color', strcolor);


strcolor = 'Black';
dim = [strx+0.28 stry-0.23 0.3 0.3];
str = '$\theta\sim 1.75^\circ$';
annotation('textbox',dim,'String',str,'LineStyle','none','FontSize',cc-5,'FontWeight', 'bold','Interpreter', 'latex','FontName','Times New Roman', 'Color', strcolor);

strcolor = 'Black';
dim = [strx+0.28 stry-0.28 0.3 0.3];
str = '\textit{Intraband}';
annotation('textbox',dim,'String',str,'LineStyle','none','FontSize',cc-5,'FontWeight', 'bold','Interpreter', 'latex','FontName','Times New Roman', 'Color', strcolor);

print('-dpdf', 'rho_vs_T_intraband.pdf');
