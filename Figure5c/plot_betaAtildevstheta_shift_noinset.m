function plot_betaAtildevstheta_shift_noinset()

figure();
hold on;
box on;

suffix = 'eh';

% %Just for legend
plot(-1 , -1, 'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 15);
plot(-1 , -1, 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 15);
plot(-1 , -1, 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 15);
plot(-1 , -1, 'ys', 'MarkerFaceColor', 'y', 'MarkerSize', 15);

%Theory
% degtheory = linspace(1.06, 2.1, 30);
% vbarmodel1 = vbarfun(degtheory);

degM = 0.8201;

degtheory = linspace(1, 2.1, 30);
% vbarmodel1 = arrayfun(@(degg) rigidvelocity(degg), degtheory);
index2 = degtheory <= 1.2295;
index2b = degtheory >= 1.23; %1.3
% vbarmodel2 = arrayfun(@(degg) relaxcarrvelocity(degg), degtheory(index2));
% vbarmodel2b = arrayfun(@(degg) relaxcarrvelocity_highangle(degg), degtheory(index2b));

index3 = degtheory <= 2;
% vbarmodel3 = arrayfun(@(degg) relaxlampsvelocity(degg), degtheory(index3));
rad = (degtheory / 180) * pi;
gamma = 1 ./ (2 * tan(rad/2));
betaA0 = 3.6;

%New
radfun = @(deg) (deg / 180) * pi;
gammafun = @(deg) 1 ./ (2 * tan(radfun(deg) / 2));
betaAtildemodel1fun = @(deg) betaA0 .* rigidvelocity(deg, degM) .* gammafun(deg);
betaAtildemodel1bfun = @(deg) betaA0 .* rigidvelocity(deg, degM);

% degtheory2 = [degtheory(index2) degtheory(index2b)];
% betaAtilde2 = [0.97 * betaAtildemodel2 smooth(smooth(betaAtildemodel2b).').'];

% betaAtildemodel2fun = @(deg) betaA0 .* relaxcarrvelocity(deg) .* gammafun(deg);
% betaAtildemodel2highdegfun = @(deg) betaA0 .* relaxcarrvelocity_highangle(deg) .* gammafun(deg);
betaAtildemodel2fun = @(deg) betaA0 .* relaxcarrvelocity_new(deg, degM) .* gammafun(deg);

betaAtildemodel3fun = @(deg) betaA0 .* relaxezzivelocity(deg, degM) .* gammafun(deg);

% %Shading
% [px, py] = plotshaded(degtheory, [betaAtildemodel1; betaAtilde2], [0.8 0.8 0.8], '');    

%Plot experiment
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


% text{1} = 'electron';
% text{2} = 'hole';

%Shading
[px, py] = plotshaded([1+4e-3 1.03 * max(degexp)-4e-3], [-1-1e-2 -1-1e-2; 1+1e-2 1+1e-2], [0.8 0.8 0.8], '');


for k = 1:2
    for deg = degexp
        if deg < 1.59 && deg ~= 1.07
            suf = suffix(k);
        elseif deg == 1.07
            suf = suffix(1); %exception
        else            
            suf = suffix(2); %exception
        end    
        load(['parameph_deg' num2str(deg) '_' suf '.mat']);
        
%         rad = (deg / 180) * pi;
%         gamma = 1 ./ (2 * tan(rad/2));
%         betaAtilde = param(1) * param(2) * gamma;
%         deltabetaAtilde = (param(1) * paramdelta(2) + paramdelta(1) * param(2)) * gamma;
        betaAtilde = param(1);
        deltabetaAtilde = paramdelta(1); 
        
%        errorbar(deg, betaAtilde, deltabetaAtilde, [clr(k) 's'], 'LineWidth', 2, 'MarkerSize', 15);
%         w = 0.04;
%         betaAtilde
%         deltabetaAtilde
%         betaAtildemodel2fun(deg)
%         if deg <= 1.2295
%             y2 = (betaAtilde + deltabetaAtilde - betaAtildemodel2fun(deg)) / betaAtildemodel2fun(deg);
%             y1 = (betaAtilde - deltabetaAtilde - betaAtildemodel2fun(deg)) / betaAtildemodel2fun(deg);
%             y = (betaAtilde - betaAtildemodel2fun(deg)) / betaAtildemodel2fun(deg);
%         else
%             y2 = (betaAtilde + deltabetaAtilde - betaAtildemodel2highdegfun(deg)) / betaAtildemodel2highdegfun(deg);
%             y1 = (betaAtilde - deltabetaAtilde - betaAtildemodel2highdegfun(deg)) / betaAtildemodel2highdegfun(deg);           
%             y = (betaAtilde - betaAtildemodel2highdegfun(deg)) / betaAtildemodel2highdegfun(deg);           
%         end

        y2 = (betaAtilde + deltabetaAtilde - betaAtildemodel2fun(deg)) / betaAtildemodel2fun(deg);
        y1 = (betaAtilde - deltabetaAtilde - betaAtildemodel2fun(deg)) / betaAtildemodel2fun(deg);
        y = (betaAtilde - betaAtildemodel2fun(deg)) / betaAtildemodel2fun(deg);
        
%         h = abs(y2 - y1);
%         w = 0.09;
%         
%         xpos = deg + [-(w / 2)  (w / 2)];
%         xposrect = (deg - w / 2);
%         wrect = w;                
        
%         k0 = 0;
%         %Exception
%         if deg == 1.06 && k == 2
%             xpos = deg + [-(w / 2)  0];
%             xposrect = deg - w / 2;
%             wrect = w / 2;
%         elseif (deg == 1.11 || deg == 1.24 || deg == 1.59 || deg == 2.02) && k == 1
%             xpos = deg + [-(w / 2)  0];
%             xposrect = deg - w / 2;
%             wrect = w / 2;
%         elseif (deg == 1.11 || deg == 1.24 || deg == 1.59 || deg == 2.02) && k == 2
%             xpos = deg + [0 (w / 2)];
%             xposrect = deg;
%             wrect = w / 2;  
%         elseif deg == 1.07 
%             xpos = deg + [-(w / 2)  (w / 2)];
%             xposrect = deg - w / 2;
%             wrect = w;     
%             k0 = k0 + 2;
%         end                
% 
%         line(xpos, y * [1 1], 'Color', clrcode{k + k0}, 'LineWidth', 4);
%         rectangle('Position', [xposrect min([y1 y2]) wrect h], 'FaceColor', [clrcode{k + k0} 0.4], 'EdgeColor', 'none', 'LineWidth', 1);    

        deltadeg = 0.018;
        if k == 1
            errorbar(deg + deltadeg, y, y2 - y, 's', 'Color', clrcode{k}, 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor', clrcode{k}); %, 'MarkerFaceColor', clrcode{k}
        else           
            errorbar(deg - deltadeg, y, y2 - y, 'o', 'Color', clrcode{k}, 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor', clrcode{k});
        end
        
        %Ezzi model
        betaAtildeth = betaAtildemodel3fun(deg);

        y2 = (betaAtilde + deltabetaAtilde - betaAtildeth) / betaAtildeth;
        y1 = (betaAtilde - deltabetaAtilde - betaAtildeth) / betaAtildeth;
        y = (betaAtilde - betaAtildeth) / betaAtildeth;
%         h = abs(y2 - y1);
% 
%         xpos = deg + [-(w / 2) (w / 2)];
%         xposrect = (deg - w / 2);
%         wrect = w;                
%                 
%         %Exception
%         k0 = 0
%         if deg == 1.06 && k == 1
%             xpos = deg + [0  (w / 2)];
%             xposrect = deg;
%             wrect = w / 2;
%         elseif (deg == 1.11 || deg == 1.24 || deg == 1.59 || deg == 2.02) && k == 1
%             xpos = deg + [-(w / 2)  0];
%             xposrect = deg - w / 2;
%             wrect = w / 2;
%         elseif (deg == 1.11 || deg == 1.24 || deg == 1.59 || deg == 2.02) && k == 2
%             xpos = deg + [0 (w / 2)];
%             xposrect = deg;
%             wrect = w / 2;      
%         elseif deg == 1.07 
%             xpos = deg + [-(w / 2)  (w / 2)];
%             xposrect = deg - w / 2;
%             wrect = w;      
%             k0 = 2;
%         end                
%         
%         line(xpos, y * [1 1], 'Color', clrcodeezzi{k + k0}, 'LineWidth', 4);
%         rectangle('Position', [xposrect min([y1 y2]) wrect h], 'FaceColor', [clrcodeezzi{k + k0} 0.4], 'EdgeColor', 'none', 'LineWidth', 1);  %clrcodeezzi{k}      

        deltadeg = 0.018;
        if k == 1
            errorbar(deg + deltadeg, y, y2 - y, 's', 'Color', clrcodeezzi{k}, 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor', clrcodeezzi{k}); %, 'MarkerFaceColor', clrcodeezzi{k}
        else
            errorbar(deg - deltadeg, y, y2 - y, 'o', 'Color', clrcodeezzi{k}, 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor', clrcodeezzi{k});
        end

    end    
end

%Plot theory

% plot(degtheory, betaAtildemodel1, 'k', 'LineWidth', 2);
% plot(degtheory2, betaAtilde2, 'b', 'LineWidth', 2);
%plot(degtheory(index3), betaAtildemodel3, 'r', 'LineWidth', 2);
% plot(degtheory, betaAtildemodel1b, 'm', 'LineWidth', 2);

xlabel('$\theta~(^\circ)$', 'FontSize', 30, 'Interpreter', 'latex');
ylabel('$(\tilde{\beta}_A^\mathrm{exp} - \tilde{\beta}_A^\mathrm{theory})/\tilde{\beta}_A^\mathrm{theory}$', 'FontSize', 30, 'Interpreter', 'latex');

% hMarkers1 = hLine1.MarkerHandle;
% hMarkers2 = hLine2.MarkerHandle;
% 
% hMarkers1.FaceColorData = uint8(255 * [0; 0; 1; 0.5]);
% hMarkers2.FaceColorData = uint8(255 * [1; 0; 0; 0.5]);
% legend('Electron - Relaxation model 1', 'Holes - Relaxation model 1', 'Electron - Relaxation model 2', 'Holes - Relaxation model 2', 'FontSize', 15, 'Interpreter', 'latex', 'Position', [.63 .67 .1 .1]);
% legend boxoff;

annotation('textbox', [.15 .8 .1 .1], 'String', {'(c) Electron-phonon theory'}, 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 30); % [.21 .8 .1 .1]
% annotation('textbox', [.16 .8 .1 .1], 'String', '(c)', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 25);
% annotation('textbox', [0.6635 0.4122 0.5304 0.0810], 'String', {'Rigid', '(geometry',' enhancement)'}, 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 20);
% annotation('textbox', [0.6510 0.8010 0.2783 0.0810], 'String', {'Relaxed model'}, 'Color', 'b', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 20);
%annotation('textbox', [0.5796 0.4 .1 .1], 'String', {'Relaxed model 2'}, 'Color', 'r', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 20);
% annotation('textbox', [0.3600 0.1762 0.5800 0.0810], 'String', {'Rigid (no geometry enhancement)'}, 'Color', 'm', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 20)

annotation('textbox', [0.59 0.23 .1 .1], 'String', 'Bad agreement', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 23); %{'Electron-phonon','theory'} [0.59 0.71 .1 .1]
annotation('textbox', [0.59 0.17 .1 .1], 'String', 'Good agreement', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 23); %{'Electron-phonon','theory'} [0.59 0.64 .1 .1]


xlim([1 1.03 * max(degexp)]);
% ylim0 = ylim;
ylim([-1 15]);

set(gca, 'LineWidth', 3, 'FontSize', 25, 'FontWeight', 'bold');

% %Inset%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axes('Position',[.46 .3 .42 .32])
% box on
% hold on
% 
% clrrigid{1} = [1 .5 0];
% clrrigid{2} = [0.34 0.7686 0.6784];
% clrrigid_nogeom{1} = [1 0 1];
% clrrigid_nogeom{2} = [.5 0 .5];
% 
% % %Just for legend
% plot(-1 , -1, 's', 'Color', clrrigid{1}, 'MarkerFaceColor', clrrigid{1}, 'MarkerSize', 15);
% plot(-1 , -1, 's', 'Color', clrrigid{2}, 'MarkerFaceColor', clrrigid{2}, 'MarkerSize', 15);
% plot(-1 , -1, 's', 'Color', clrrigid_nogeom{1}, 'MarkerFaceColor', clrrigid_nogeom{1}, 'MarkerSize', 15);
% plot(-1 , -1, 's', 'Color', clrrigid_nogeom{2}, 'MarkerFaceColor', clrrigid_nogeom{2}, 'MarkerSize', 15);
% 
% for k = 1:2
%     for deg = degexp
%         if deg < 1.59
%             suf = suffix(k);
%         else
%             suf = suffix(2); %exception
%         end    
%         load(['parameph_deg' num2str(deg) '_' suf '.mat']);
%         
%         betaAtilde = param(1);
%         deltabetaAtilde = paramdelta(1); 
%         
%         betaAtildeth = betaAtildemodel1fun(deg);
%         y2 = (betaAtilde + deltabetaAtilde - betaAtildeth) / betaAtildeth;
%         y1 = (betaAtilde - deltabetaAtilde - betaAtildeth) / betaAtildeth;
%         y = (betaAtilde - betaAtildeth) / betaAtildeth;
%         
%         h = abs(y2 - y1);
%         w = 0.06;
%         [(deg - w / 2) min([y1 y2]) w h]
% 
%         line(deg + [-(w / 2)  (w / 2)], y * [1 1], 'Color', clrrigid{k}, 'LineWidth', 2);
%         rectangle('Position', [(deg - w / 2) min([y1 y2]) w h], 'FaceColor', [clrrigid{k} 0.5], 'EdgeColor', 'none', 'LineWidth', 1);            
%         
%         %No geometry enhancement
%         betaAtildeth = betaAtildemodel1bfun(deg);
%         y2 = (betaAtilde + deltabetaAtilde - betaAtildeth) / betaAtildeth;
%         y1 = (betaAtilde - deltabetaAtilde - betaAtildeth) / betaAtildeth;
%         y = (betaAtilde - betaAtildeth) / betaAtildeth;
%         
%         h = abs(y2 - y1);
%         [(deg - w / 2) min([y1 y2]) w h]
% 
%         line(deg + [-(w / 2)  (w / 2)], y * [1 1], 'Color', clrrigid_nogeom{k}, 'LineWidth', 2);
%         rectangle('Position', [(deg - w / 2) min([y1 y2]) w h], 'FaceColor', [clrrigid_nogeom{k} 0.5], 'EdgeColor', 'none', 'LineWidth', 1);                    
%     end    
% end
% 
% xlabel('$\theta~(^\circ)$', 'FontSize', 15, 'Interpreter', 'latex');
% ylabel('$(\tilde{\beta}_A^\mathrm{exp} - \tilde{\beta}_A^\mathrm{theory})/\tilde{\beta}_A^\mathrm{theory}$', 'FontSize', 15, 'Interpreter', 'latex');
% 
% % xlim([1.55 2.05]);
% %ylim([0 7]); %[1.6 7]
% 
% legend('Electron - Rigid', 'Holes - Rigid', 'Electron - Rigid (no geometry enhanced)', 'Holes - Rigid (no geometry enhanced)', 'FontSize', 10, 'Interpreter', 'latex', 'Position', [.64 .49 .1 .1]);
% legend boxoff;
% 
% set(gca, 'LineWidth', 2, 'FontSize', 15, 'FontWeight', 'bold');
% 
% xlim([1 1.025 * max(degexp)]);
% ylim0 = ylim;
% ylim([min(ylim0) 1.1 * max(ylim0)]);

print('-dpdf', 'betaAtilde_vs_theta_eph_shift_noinset.pdf');

end

%Shading
    function [px, py] = plotshaded(x, y, fstr, dispname)

	px = [x,fliplr(x)]; % make closed patch
%    px = [x(1, :), fliplr(x(2, :))]; % make closed patch
	py = [y(1, :), fliplr(y(2, :))];
	p(1) = patch(px, py, 1,'FaceColor', fstr, 'EdgeColor', 'none', 'DisplayName', dispname); %'FaceAlpha', 0.2
    end