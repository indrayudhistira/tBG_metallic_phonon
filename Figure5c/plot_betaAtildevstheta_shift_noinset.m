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

degM = 0.7210; %0.8201

degtheory = linspace(1, 2.1, 30);
rad = (degtheory / 180) * pi;
% gamma = 1 ./ (2 * tan(rad/2));
betaA0 = 3.6;

%New
radfun = @(deg) (deg / 180) * pi;
gammafun = @(deg) 1 ./ (2 * tan(radfun(deg) / 2));
betaAtildemodel2fun = @(deg) betaA0 .* relaxcarrvelocity_new(deg, degM) .* gammafun(deg);
betaAtildemodel3fun = @(deg) betaA0 .* relaxezzivelocity(deg, degM) .* gammafun(deg);

%Plot experiment
degexp = [1.06 1.07 1.08 1.11 1.24 1.59 2.02];

clrcode{1} = [1 0 0];
clrcode{2} = [1 0 0];

clrcode{3} = [0 0 0];
clrcode{4} = [0.4 0.2 0.6];

clrcodeezzi{1} = [0 0 1];
clrcodeezzi{2} = [0 0 1];

clrcodeezzi{3} = [	000 106 110]/255;
clrcodeezzi{4} = [0 1 1];

%Shading
[px, py] = plotshaded([1+4e-3 1.03 * max(degexp)-4e-3], [-1-1e-2 -1-1e-2; 1+1e-2 1+1e-2], [0.8 0.8 0.8], '');

radfun = @(deg) (deg / 180) * pi;
gammafun = @(deg) 1 ./ (2 * tan(radfun(deg) / 2));

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
        
        betaAtilde = param(1);
        deltabetaAtilde = paramdelta(1); 
        

        y2 = (betaAtilde + deltabetaAtilde - betaAtildemodel2fun(deg)) / betaAtildemodel2fun(deg);
%         y1 = (betaAtilde - deltabetaAtilde - betaAtildemodel2fun(deg)) / betaAtildemodel2fun(deg);
        y = (betaAtilde - betaAtildemodel2fun(deg)) / betaAtildemodel2fun(deg);
        
        betacarr = betaAtilde / (relaxcarrvelocity_new(deg, degM) * gammafun(deg))                

        deltadeg = 0.018;
        if k == 1
            errorbar(deg + deltadeg, y, y2 - y, 's', 'Color', clrcode{k}, 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor', clrcode{k}); %, 'MarkerFaceColor', clrcode{k}
        else           
            errorbar(deg - deltadeg, y, y2 - y, 'o', 'Color', clrcode{k}, 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor', clrcode{k});
        end
        
        %Ezzi model
        betaAtildeth = betaAtildemodel3fun(deg);

        y2 = (betaAtilde + deltabetaAtilde - betaAtildeth) / betaAtildeth;
%         y1 = (betaAtilde - deltabetaAtilde - betaAtildeth) / betaAtildeth;
        y = (betaAtilde - betaAtildeth) / betaAtildeth;
        
        betaezzi = betaAtilde / (relaxezzivelocity(deg, degM) * gammafun(deg))        
        
        deltadeg = 0.018;
        if k == 1
            errorbar(deg + deltadeg, y, y2 - y, 's', 'Color', clrcodeezzi{k}, 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor', clrcodeezzi{k}); %, 'MarkerFaceColor', clrcodeezzi{k}
        else
            errorbar(deg - deltadeg, y, y2 - y, 'o', 'Color', clrcodeezzi{k}, 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor', clrcodeezzi{k});
        end

    end    
end

xlabel('$\theta~(^\circ)$', 'FontSize', 30, 'Interpreter', 'latex');
ylabel('$(\tilde{\beta}_A^\mathrm{exp} - \tilde{\beta}_A^\mathrm{theory})/\tilde{\beta}_A^\mathrm{theory}$', 'FontSize', 30, 'Interpreter', 'latex');

annotation('textbox', [.15 .8 .1 .1], 'String', {'(c) Electron-phonon theory'}, 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 30); % [.21 .8 .1 .1]

annotation('textbox', [0.59 0.23 .1 .1], 'String', 'Bad agreement', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 23); %{'Electron-phonon','theory'} [0.59 0.71 .1 .1]
annotation('textbox', [0.59 0.17 .1 .1], 'String', 'Good agreement', 'LineStyle', 'none','Interpreter', 'latex', 'FontSize', 23); %{'Electron-phonon','theory'} [0.59 0.64 .1 .1]


xlim([1 1.03 * max(degexp)]);
ylim([-1 15]);

set(gca, 'LineWidth', 3, 'FontSize', 25, 'FontWeight', 'bold');

print('-dpdf', 'Figure5c.pdf');

end

%Shading
    function [px, py] = plotshaded(x, y, fstr, dispname)

	px = [x,fliplr(x)]; % make closed patch
%    px = [x(1, :), fliplr(x(2, :))]; % make closed patch
	py = [y(1, :), fliplr(y(2, :))];
	p(1) = patch(px, py, 1,'FaceColor', fstr, 'EdgeColor', 'none', 'DisplayName', dispname); %'FaceAlpha', 0.2
    end