function plot_vFvstheta_shift_noinset()

figure();
hold on;
box on;

suffix = 'eh';
clr = 'br';

% %Just for legend
hLine1 = plot(-1 , -1, 'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 15);
hLine2 = plot(-1 , -1, 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 15);
hLine3 = plot(-1 , -1, 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 15);
hLine4 = plot(-1 , -1, 'ys', 'MarkerFaceColor', 'y', 'MarkerSize', 15);

degexp = [1.06 1.07 1.08 1.11 1.24 1.59 2.02];

clrcode{1} = [1 0 0];
% clrcode{2} = [1 0 0];
clrcode{2} = [1 0 0];

clrcode{3} = [0 0 0];
clrcode{4} = [0.4 0.2 0.6];

clrcodeezzi{1} = [0 0 1];
% clrcodeezzi{2} = [1 1 0];
clrcodeezzi{2} = [0 0 1];

clrcodeezzi{3} = [	000 106 110]/255;
clrcodeezzi{4} = [0 1 1];

degM = 0.8201;

%Shading
[px, py] = plotshaded([1+4e-3 1.03 * max(degexp)-4e-3], [-1-1e-2 -1-1e-2; 1+1e-2 1+1e-2], [0.8 0.8 0.8], '');

%Plot experiment
for k = 1:2
    for deg = degexp
        if deg < 1.59 && deg ~= 1.07
            suf = suffix(k);
        elseif deg == 1.07
            suf = suffix(1); %exception
        else            
            suf = suffix(1); %exception
        end   
        load(['paramPlanckian_deg' num2str(deg) '_' suf '.mat']);

%         if deg <= 1.2295
%             y2 = (param(2) + paramdelta(2) - relaxcarrvelocity(deg)) / relaxcarrvelocity(deg);
%             y1 = (param(2) - paramdelta(2) - relaxcarrvelocity(deg)) / relaxcarrvelocity(deg);
%             y = (param(2) - relaxcarrvelocity(deg)) / relaxcarrvelocity(deg);
%         else
%             y2 = (param(2) + paramdelta(2) - relaxcarrvelocity_highangle(deg)) / relaxcarrvelocity_highangle(deg);
%             y1 = (param(2) - paramdelta(2) - relaxcarrvelocity_highangle(deg)) / relaxcarrvelocity_highangle(deg);           
%             y = (param(2) - relaxcarrvelocity_highangle(deg)) / relaxcarrvelocity_highangle(deg);
%         end
        
        vFth = relaxcarrvelocity_new(deg, degM);
        y2 = (param(2) + paramdelta(2) - vFth) / vFth;
        y1 = (param(2) - paramdelta(2) - vFth) / vFth;           
        y = (param(2) - vFth) / vFth;
                
%         h = abs(y2 - y1);
%         w = 0.09;
        
%         ['paramPlanckian_deg' num2str(deg) '_' suf '.mat']
%         deg
%         param(2)
%         param
%         vFth        

%         xpos = deg + [-(w / 2)  (w / 2)];
%         xposrect = (deg - w / 2);
%         wrect = w;
% 
%         %Exception
%         k0 = 0;
%         if (deg == 1.06 && k == 2)
%             xpos = deg + [-(w / 2)  0];
%             wrect = w / 2;
%         elseif (deg == 1.11 || deg == 1.24 || deg == 1.59 || deg == 2.02) && k == 1
%             xpos = deg + [-(w / 2)  -(w / 4)];            
%             wrect = w / 4;            
%         elseif (deg == 1.11 || deg == 1.24 || deg == 1.59 || deg == 2.02) && k == 2
%             xpos = deg + [-(w / 4) 0];
%             xposrect = deg - w / 4;
%             wrect = w / 4;
%         elseif deg == 1.07
%             xpos = deg + [-(w / 2) (w / 2)];
%             xposrect = deg - w / 2;
%             wrect = w;            
%             k0 = 2;
%         end
%         
%         line(xpos, y * [1 1], 'Color', clrcode{k + k0}, 'LineWidth', 4);
%         rectangle('Position', [xposrect min([y1 y2]) wrect h], 'FaceColor', [clrcode{k + k0} 0.4], 'LineWidth', 1, 'EdgeColor', 'none'); %, 'EdgeColor', clrcode{k}
        
%         errorbar(deg, param(2), paramdelta(2), [clr(k) 's'], 'LineWidth', 2, 'MarkerSize', 10);

        deltadeg = 0.018;
        if k == 1
            errorbar(deg + deltadeg, y, y2 - y, 's', 'Color', clrcode{k}, 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor', clrcode{k}); %
        else           
            errorbar(deg - deltadeg, y, y2 - y, 'o', 'Color', clrcode{k}, 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor', clrcode{k});
        end
        
        %Ezzi model
        vFtheory = relaxezzivelocity(deg, degM);
        y2 = (param(2) + paramdelta(2) - vFtheory) / vFtheory;
        y1 = (param(2) - paramdelta(2) - vFtheory) / vFtheory;
        y = (param(2) - vFtheory) / vFtheory;
%         h = abs(y2 - y1);
% 
%         xpos = deg + [-(w / 2)  (w / 2)];
%         xposrect = (deg - w / 2);
%         wrect = w;        
%         
%         %Exception
%         k0 = 0;
%         if deg == 1.06 && k == 1
%             xpos = deg + [0  (w / 2)];
%             xposrect = deg;
%             wrect = w / 2;
%         elseif (deg == 1.11 || deg == 1.24 || deg == 1.59 || deg == 2.02) && k == 1
%             xpos = deg + [0 (w / 4)];
%             xposrect = deg;
%             wrect = w / 4;
%         elseif (deg == 1.11 || deg == 1.24 || deg == 1.59 || deg == 2.02) && k == 2
%             xpos = deg + w / 4 + [0 (w / 4)];
%             xposrect = deg + w / 4;
%             wrect = w / 4;          
%         elseif deg == 1.07
%             xpos = deg + [-(w / 2) (w / 2)];
%             xposrect = deg - w / 2;
%             wrect = w;            
%             k0 = 2;            
%         end        
%         
%         line(xpos, y * [1 1], 'Color', clrcodeezzi{k + k0}, 'LineWidth', 4);
%         rectangle('Position', [xposrect min([y1 y2]) wrect h], 'FaceColor', [clrcodeezzi{k + k0} 0.4], 'EdgeColor', 'none', 'LineWidth', 1);  %clrcodeezzi{k}      

        deltadeg = 0.018;
        if k == 1
            errorbar(deg + deltadeg, y, y2 - y , 's', 'Color', clrcodeezzi{k}, 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor', clrcodeezzi{k}); %
        else           
            errorbar(deg - deltadeg, y, y2 - y, 'o', 'Color', clrcodeezzi{k}, 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor', clrcodeezzi{k});
        end

    end    
end

%Plot theory
% degtheory = linspace(1.06, 2.1, 30);
% vbarmodel1 = vbarfun(degtheory);
% plot(degtheory, vbarmodel1, 'k', 'LineWidth', 2);
% index2 = degtheory <= 1.2295;
% index2b = degtheory >= 1.3;
% vbarmodel2 = arrayfun(@(degg) relaxcarrvelocity(degg), degtheory(index2));
% vbarmodel2b = arrayfun(@(degg) relaxcarrvelocity_highangle(degg), degtheory(index2b));
% index3 = degtheory <= 2;
% vbarmodel3 = arrayfun(@(degg) relaxlampsvelocity(degg), degtheory(index3));
% plot(degtheory, vbarmodel1, 'k', 'LineWidth', 2);
% plot([degtheory(index2) degtheory(index2b)], [vbarmodel2 smooth(smooth(vbarmodel2b).').'], 'b', 'LineWidth', 2);
% plot(degtheory(index3), vbarmodel3, 'r', 'LineWidth', 2);

%Theory
% degtheory = linspace(1.06, 2.1, 30);
% vbarmodel1 = vbarfun(degtheory);
degtheory = linspace(1, 2.1, 30);
% vbarmodel1 = arrayfun(@(degg) rigidvelocity(degg), degtheory);
% plot(degtheory, vbarmodel1, 'k', 'LineWidth', 2);
index2 = degtheory <= 1.2295;
index2b = degtheory >= 1.23; %1.3
% vbarmodel2 = arrayfun(@(degg) relaxcarrvelocity(degg), degtheory(index2));
% vbarmodel2b = arrayfun(@(degg) relaxcarrvelocity_highangle(degg), degtheory(index2b));
% index3 = degtheory <= 2;
% vbarmodel3 = arrayfun(@(degg) relaxlampsvelocity(degg), degtheory(index3));

degtheory2 = [degtheory(index2) degtheory(index2b)];
% vbar2 = [0.968*vbarmodel2 smooth(smooth(vbarmodel2b).').'];

% %Shading
% [px, py] = plotshaded(degtheory, [vbarmodel1; vbar2], [0.8 0.8 0.8], '');

%Plot theory
% plot(degtheory, vbarmodel1, 'k', 'LineWidth', 2);
% plot(degtheory2, vbar2, 'b', 'LineWidth', 2);
%plot(degtheory(index3), vbarmodel3, 'r', 'LineWidth', 2);

xlabel('$\theta~(^\circ)$', 'FontSize', 30, 'Interpreter', 'latex');
% ylabel('$\Delta v_F/v_F^\mathrm{theory}$', 'FontSize', 30, 'Interpreter', 'latex');
ylabel('$(v_F^\mathrm{exp} - v_F^\mathrm{theory})/v_F^\mathrm{theory}$', 'FontSize', 30, 'Interpreter', 'latex');

% legend('Electron - Relaxation model 1', 'Holes - Relaxation model 1', 'Electron - Relaxation model 2', 'Holes - Relaxation model 2', 'FontSize', 15, 'Interpreter', 'latex', 'Position', [.29 .68 .1 .1]);
% legend boxoff;

annotation('textbox', [0.14 0.80 .1 .1], 'String', '(b)', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 30);
annotation('textbox', [0.23 0.79 .1 .1], 'String', {'Planckian', 'theory'}, 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 30);
% annotation('textbox', [.15 .82 .1 .1], 'String', '', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 25);
% a=annotation('textbox', [0.6710 0.48 .1 .1], 'String', {'Rigid lattice'}, 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 22);
% b=annotation('textbox', [0.5532 0.62 .1 .1], 'String', {'Relaxation model 1'}, 'Color', 'b', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 22);
% c=annotation('textbox', [0.5461 0.24 .1 .1], 'String', {'Relaxation model 2'}, 'Color', 'r', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 22);

annotation('textbox', [0.56 0.26 .1 .1], 'String', 'Bad agreement', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 23); %{'Electron-phonon','theory'}
annotation('textbox', [0.56 0.195 .1 .1], 'String', 'Good agreement', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 23); %{'Electron-phonon','theory'}


xlim([1 1.03 * max(degexp)]);
ylim0 = ylim;
% ylim([min(ylim0) 1.1 * max(ylim0)]);
ylim([-1 11]);
% ylim([01.7]);

set(gca, 'LineWidth', 3, 'FontSize', 25, 'FontWeight', 'bold');

% %Inset%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axes('Position',[.32 .31 .29 .32]);
% box on
% hold on
% 
% clrrigid{1} = [1 .5 0];
% clrrigid{2} = [0.34 0.7686 0.6784];
% 
% %Just for legend
% plot(-1 , -1, 's', 'Color', clrrigid{1}, 'MarkerFaceColor', clrrigid{1}, 'MarkerSize', 15);
% plot(-1 , -1, 's', 'Color', clrrigid{2}, 'MarkerFaceColor', clrrigid{2}, 'MarkerSize', 15);
% 
% for k = 1:2
%     for deg = degexp
%         if deg < 1.59
%             suf = suffix(k);
%         else
%             suf = suffix(1); %exception
%         end
%         
%         load(['paramPlanckian_deg' num2str(deg) '_' suf '.mat']);
%         
%         vFth = rigidvelocity(deg);
%         y2 = (param(2) + paramdelta(2) - vFth) / vFth;
%         y1 = (param(2) - paramdelta(2) - vFth) / vFth;
%         y = (param(2) - vFth) / vFth;
%         
%         ['paramPlanckian_deg' num2str(deg) '_' suf '.mat']
%         deg
%         param(2)
%         param
%         vFth
%         
%         h = abs(y2 - y1);
%         w = 0.06;
%         [(deg - w / 2) min([y1 y2]) w h]
% 
%         line(deg + [-(w / 2)  (w / 2)], y * [1 1], 'Color', clrrigid{k}, 'LineWidth', 2);
%         rectangle('Position', [(deg - w / 2) min([y1 y2]) w h], 'FaceColor', [clrrigid{k} 0.5], 'EdgeColor', 'none', 'LineWidth', 1);                    
%     end    
% end
% 
% xlabel('$\theta~(^\circ)$', 'FontSize', 15, 'Interpreter', 'latex');
% ylabel('$(v_F^\mathrm{exp} - v_F^\mathrm{theory})/v_F^\mathrm{theory}$', 'FontSize', 15, 'Interpreter', 'latex');
% 
% legend('Electron - Rigid', 'Holes - Rigid', 'FontSize', 13, 'Interpreter', 'latex', 'Position', [.3 .53 .1 .1]);
% legend boxoff;
% 
% set(gca, 'LineWidth', 2, 'FontSize', 15, 'FontWeight', 'bold');
% 
% xlim([1 1.02 * max(degexp)]);
% ylim0 = ylim;
% % ylim([.2 * min(ylim0) 1.4 * max(ylim0)]);

print('-dpdf', 'vF_vs_theta_Planckian_shift_noinset.pdf');

end

%Shading
    function [px, py] = plotshaded(x, y, fstr, dispname)

	px = [x,fliplr(x)]; % make closed patch
%    px = [x(1, :), fliplr(x(2, :))]; % make closed patch
	py = [y(1, :), fliplr(y(2, :))];
	p(1) = patch(px, py, 1,'FaceColor', fstr, 'EdgeColor', 'none', 'DisplayName', dispname); %'FaceAlpha', 0.2
    end