function plot_Cvstheta()

%clear;

figure();
hold on;
box on;

suffix = 'eh';

clrcode{1} = [1 .5 0.1725];
clrcode{2} = [1 0 1];

hLine1 = plot(-1 , -2, 's', 'Color', clrcode{1}, 'MarkerSize', 15, 'LineWidth', 2, 'MarkerFaceColor', clrcode{1}); %'MarkerFaceColor', clrcode{1}
hLine2 = plot(-1 , -2, 'o', 'Color', clrcode{2}, 'MarkerSize', 15, 'LineWidth', 2, 'MarkerFaceColor', clrcode{2}); %, 'MarkerFaceColor', clrcode{2}

degexp = [1.06 1.07 1.08 1.11 1.24 1.59 2.02];

%Shading
[px, py] = plotshaded([1+4e-3 1.03 * max(degexp)-4e-3], [1-1e-2 1-1e-2; 0+1e-2 0+1e-2], [0.8 0.8 0.8], '');    

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
        
        y = param(1);
        y2 = y + paramdelta(1);
        y1 = y - paramdelta(1);
      
        deltadeg = 0.018;
        if k == 1
            errorbar(deg + deltadeg, y, y2 - y, 's', 'Color', clrcode{k}, 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor', clrcode{k}); %, 'MarkerFaceColor', clrcode{k}
        else           
            errorbar(deg - deltadeg, y, y2 - y, 'o', 'Color', clrcode{k}, 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor', clrcode{k});
        end
    end    
end

%Plot theory
degtheory = linspace(1, 2.1, 30);
Cminonetheory = zeros(size(degtheory));
plot(degtheory, Cminonetheory, 'k', 'LineWidth', 2);

xlabel('$\theta~(^\circ)$', 'FontSize', 30, 'Interpreter', 'latex');
ylabel('$C$', 'FontSize', 30, 'Interpreter', 'latex');

legend('Electron', 'Holes', 'FontSize', 18, 'Interpreter', 'latex'); %, 'Position', [.6 .67 .1 .1]
% legend boxoff;

annotation('textbox', [.14 .8 .1 .1], 'String', '(d) Planckian theory', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 30);
annotation('textbox', [0.58 0.185 .1 .1], 'String', 'Bad agreement', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 23); %{'Electron-phonon','theory'} [0.59 0.39 .1 .1]
annotation('textbox', [0.48 0.135 .1 .1], 'String', 'Good agreement', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 23); %{'Electron-phonon','theory'} [0.59 0.33 .1 .1]

xlim([1 1.03 * max(degexp)]);
ylim([0 15]);

set(gca, 'LineWidth', 3, 'FontSize', 25, 'FontWeight', 'bold');

print('-dpdf', 'C_vs_theta_Planckian.pdf');

end

%Shading
    function [px, py] = plotshaded(x, y, fstr, dispname)

	px = [x,fliplr(x)]; % make closed patch
	py = [y(1, :), fliplr(y(2, :))];
	p(1) = patch(px, py, 1,'FaceColor', fstr, 'EdgeColor', 'none', 'DisplayName', dispname); %'FaceAlpha', 0.2
    end