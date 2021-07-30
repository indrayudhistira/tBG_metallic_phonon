%%
clc
clear all

load Rho_vs_theta_rigid_new_T10_n5.mat
cc=25;
figure;


rhoephns = abs(rhoephns);
[ind1, ind2] = min(rhoephns);

%Plot vertical line where va = vf
%ymax = 0.07; %0.042 0.25
%plot([deg_crit deg_crit], [0 ymax], 'k--', 'LineWidth', 2);

%Plot vertical line where theta = theta_M
%plot([deg_magic deg_magic], [0 ymax], 'k--', 'LineWidth', 2);


%yyaxis left
% 145 is position of magic angle in relaxcarr
% 91 is position of magic angle in rigid
pl = plot(degg, rhoephns ,'r-', 'LineWidth', 3);
hold on ;
%pl2 = plot(deg(1:cutdegval), log(normresrigidblochgrun(1:cutdegval)));

%yyaxis right
%pr=plot(deg,Data_v_relax(:,2),'b','LineWidth',3);

% upp = 5;
ycross = linspace(0,.5,1000);
% ycross2 = 10.^linspace(0,6,1000);
for ii = 1:length(cross)
    plot(linspace(cross(ii),cross(ii),1000),ycross,'-.k','LineWidth',2)
    
end

plot(linspace(cross(1),cross(1),10),linspace(0.003,0.01,10),'r','LineWidth',3)
plot(linspace(cross(2),cross(2),10),linspace(0.003,0.01,10),'r','LineWidth',3)

thetam = 1.05;
plot(linspace(thetam,thetam,1000),ycross,'-.b','LineWidth',2)

%cphline = linspace(1.623e4/0.8739e6,1.623e4/0.8739e6,91);
%plot(deg,cphline,'k')
%xlabel('$n (10^{10} \mathrm{cm}^{-2} )$', 'FontSize', 30, 'Interpreter', 'latex');
xlabel('$\theta (^\circ)$', 'FontSize', cc, 'Interpreter', 'latex');
%yyaxis left
ylabel('$\rho_{e-ph} (h/e^2)$', 'FontSize', cc, 'Interpreter', 'latex');
%yyaxis right
%ylabel('$v_F/v_0$', 'FontSize', 30, 'Interpreter', 'latex','Color','blue');

%l = legend(p);
%set(l, 'FontSize', 18,'Location', 'NorthWest', 'Interpreter', 'latex');
%legend boxoff;

txstr{1} = ['Rigid model'];
txstr{2} = ['$n=' sprintf('%g', n) '\times 10^{10}~\mathrm{cm}^{-2}$'];
txstr{3} = ['$T=' sprintf('%g', T)  '$~K'];
annotation('textbox', [0.25 0.8 0.1 0.1], 'String', txstr, 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', cc-2);

% text(1.04*deg_magic, 0.01, '$\theta_M$','Interpreter', 'latex', 'FontSize', 22);
% text(1.05*deg_crit, 0.007, '$\theta_{cr}$','Interpreter', 'latex', 'FontSize', 22);


xx1 = linspace(0.3, cross(1),10);
yy1 = linspace(0,5,10);
a11=area(xx1,yy1);
a11.FaceAlpha = 0.1;
a11.FaceColor = [1 0 0];
%
xx1 = linspace(cross(1), cross(2),10);
yy1 = linspace(0,10,10);
a11=area(xx1,yy1);
a11.FaceAlpha = 0.19;
a11.FaceColor = [0 0 1];

xx1 = linspace(cross(2), cross(3),10);
yy1 = linspace(0,100,10);
a11=area(xx1,yy1);
a11.FaceAlpha = 0.1;
a11.FaceColor = [1 0 0];
%
xx1 = linspace(cross(3), cross(4),10);
yy1 = linspace(0,100,10);
a11=area(xx1,yy1);
a11.FaceAlpha = 0.19;
a11.FaceColor = [0 0 1];


xx1 = linspace(cross(4), 1.21,10);
yy1 = linspace(0,20,10);
a11=area(xx1,yy1);
a11.FaceAlpha = 0.1;
a11.FaceColor = [1 0 0];



h = gca;
h.XMinorTick='on';
h.YMinorTick='on';

xlim([.42 1.21]);
ylim([-0 0.31]);

stry=0.73;
dim = [0.36 stry-0.45 0.3 0.3];
str = {'$v_F>c_{ph}$' };
annotation('textbox',dim,'String',str,'LineStyle','none','FontSize',cc-4,'FontWeight', 'bold','Interpreter', 'latex','FontName','Times New Roman', 'Color', 'Black');


dim = [0.58 stry-0.45 0.3 0.3];
str = {'$v_F<c_{ph}$' };
annotation('textbox',dim,'String',str,'LineStyle','none','FontSize',cc-4,'FontWeight', 'bold','Interpreter', 'latex','FontName','Times New Roman', 'Color', 'Black');


dim = [0.52 stry-0.73 0.3 0.3];
str = {'$\theta_{cr}$' };
annotation('textbox',dim,'String',str,'LineStyle','none','FontSize',cc-4,'FontWeight', 'bold','Interpreter', 'latex','FontName','Times New Roman', 'Color', 'Black');

dim = [0.7 stry-0.73 0.3 0.3];
str = {'$\theta_M$' };
annotation('textbox',dim,'String',str,'LineStyle','none','FontSize',cc-4,'FontWeight', 'bold','Interpreter', 'latex','FontName','Times New Roman', 'Color', 'Black');


dim = [0.81 stry-0.73 0.3 0.3];
str = {'$\theta_{cr}$' };
annotation('textbox',dim,'String',str,'LineStyle','none','FontSize',cc-4,'FontWeight', 'bold','Interpreter', 'latex','FontName','Times New Roman', 'Color', 'Black');



set(gca, 'LineWidth', 3, 'FontSize', cc, 'FontWeight', 'bold');

% yticks([0 .1 .2 .3]);
% yticklabels({'0', '0.1',  '0.2', '0.3'});

annotation('textbox', [0 0.8 .1 .1], 'String', '(a)', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 25);

box on;


ax = gca;
%exportgraphics(ax,'Fig-Rigid.pdf','Resolution',300) 

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
    print('-dpdf', 'Fig-Rigid.pdf');

    % Restore the previous settings
    set(gcf,'PaperType',prePaperType);
    set(gcf,'PaperUnits',prePaperUnits);
    set(gcf,'Units',preUnits);
    set(gcf,'PaperPosition',prePaperPosition);
    set(gcf,'PaperSize',prePaperSize);    


