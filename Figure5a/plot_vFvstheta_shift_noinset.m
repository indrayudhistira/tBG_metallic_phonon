function plot_vFvstheta_shift_noinset()

figure();
hold on;
box on;

suffix = 'eh';
% clr = 'br';

hLine1 = plot(-1 , -1, 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 15);
hLine2 = plot(-1 , -1, 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 15);
hLine3 = plot(-1 , -1, 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 15);
hLine4 = plot(-1 , -1, 'ys', 'MarkerFaceColor', 'y', 'MarkerSize', 15);

degexp = [1.06 1.07 1.08 1.11 1.24 1.59 2.02];

%Theory
degM = 0.7210; %0.8201

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


text{1} = 'electron';
text{2} = 'hole';

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
            suf = suffix(2); %exception
        end       
        load(['parameph_deg' num2str(deg) '_' suf '.mat']);       

        vFth = relaxcarrvelocity_new(deg, degM);
        y2 = (param(2) + paramdelta(2) - vFth) / vFth;
%         y1 = (param(2) - paramdelta(2) - vFth) / vFth;
        y = (param(2) - vFth) / vFth;

        deltadeg = 0.018;
        if k == 1
            errorbar(deg + deltadeg, y, y2 - y, 's', 'Color', clrcode{k}, 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor', clrcode{k}); %, 'MarkerFaceColor', clrcode{k}
        else           
            errorbar(deg - deltadeg, y, y2 - y, 'o', 'Color', clrcode{k}, 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor', clrcode{k});
        end
        
        %Ezzi model
        vFth = relaxezzivelocity(deg, degM);

        y2 = (param(2) + paramdelta(2) - vFth) / vFth;
%         y1 = (param(2) - paramdelta(2) - vFth) / vFth;
        y = (param(2) - vFth) / vFth;
%         h = abs(y2 - y1);
 
        deltadeg = 0.018;
        if k == 1
            errorbar(deg + deltadeg, y, y2 - y, 's', 'Color', clrcodeezzi{k}, 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor', clrcodeezzi{k}); %, 'MarkerFaceColor', clrcodeezzi{k}
        else
            errorbar(deg - deltadeg, y, y2 - y, 'o', 'Color', clrcodeezzi{k}, 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor', clrcodeezzi{k});
        end
    end    
end

xlabel('$\theta~(^\circ)$', 'FontSize', 30, 'Interpreter', 'latex');
% ylabel('$v_F/v_0$', 'FontSize', 30, 'Interpreter', 'latex');
ylabel('$(v_F^\mathrm{exp} - v_F^\mathrm{theory})/v_F^\mathrm{theory}$', 'FontSize', 30, 'Interpreter', 'latex');

annotation('textbox', [0.15 0.8 .1 .1], 'String', '(a) Electron-phonon theory', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 30); %{'Electron-phonon','theory'}

annotation('textbox', [0.58 0.275 .1 .1], 'String', 'Bad agreement', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 23); %{'Electron-phonon','theory'}
annotation('textbox', [0.58 0.21 .1 .1], 'String', 'Good agreement', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 23); %{'Electron-phonon','theory'}

xlim([1 1.03 * max(degexp)]);
% ylim0 = ylim;
ylim([-1 10]);

set(gca, 'LineWidth', 3, 'FontSize', 25, 'FontWeight', 'bold');

print('-dpdf', 'Figure5a.pdf');

end

%Shading
    function [px, py] = plotshaded(x, y, fstr, dispname)

	px = [x,fliplr(x)]; % make closed patch
%    px = [x(1, :), fliplr(x(2, :))]; % make closed patch 
	py = [y(1, :), fliplr(y(2, :))];
	p(1) = patch(px, py, 1,'FaceColor', fstr, 'EdgeColor', 'none', 'DisplayName', dispname); %'FaceAlpha', 0.2
    end