clc
clear all
close all

% --------------input parameters
L_half = 30;
Delta_x = 0.25;
Pe = 4000;
InjectionPoint_normalized = 0.2;
NumRandomWalkers = 200000;
NumSteps = 100;
Factor_related_to_delatT = 10;

%------------------------------
Lx = [-L_half:Delta_x:L_half]';
LengthAll = max(Lx) - min(Lx);
InjectionPoint = LengthAll * InjectionPoint_normalized + min(Lx);

NumPnts = size(Lx, 1);
Element = [[1:size(Lx, 1) - 1]', [2:size(Lx, 1)]'];
NumSubPlot = 2;

% Velocity = 1e-10 .* ((Lx(Element(:, 1)) + Lx(Element(:, 2))) .* 0.5) .^ 4;
Velocity = rand(NumPnts-1, 1) .* 1e5;

meanV = mean(Velocity);
t_ = 1 / meanV;
deltaT = t_ / Factor_related_to_delatT;

Dm = meanV * L_half * 2 / Pe;

% ----------------for plot
Height_half = 0.5;
Pnt_y = [[Lx, zeros(NumPnts, 1) + Height_half]; [Lx, zeros(NumPnts, 1) - Height_half]];
Element_y = [[1:size(Lx, 1) - 1]', [2:size(Lx, 1)]', [NumPnts + 2:NumPnts * 2]', [NumPnts + 1:NumPnts * 2 - 1]'];

figure(1);
subplot(NumSubPlot, 1, 1)
title('Velocity distribution');
xlabel('x (m)'); ylabel('y (m)');
hold on
patch('Vertices', Pnt_y, 'Faces', Element_y, 'FaceVertexCData', Velocity, 'FaceColor', 'flat', 'EdgeAlpha', 0.0, 'facealpha', 1); hold on
xlim([min(Lx) - 2, max(Lx) + 2])
ylim([-2, 2])
pbaspect([(max(Lx) + 2) * 1, 4, 1]); hold on
hjk = colorbar;
hjk.Label.String = 'Velocity [L/T]';
%-----------------------

% ParticleID, ElementID, Position, Velocity
PointMatrix = zeros(NumRandomWalkers, 4);
PointMatrix([1:NumRandomWalkers], 1) = [1:NumRandomWalkers]';

ElementIDinit = find(Lx >= InjectionPoint, 1);
PointMatrix([1:NumRandomWalkers], 2) = PointMatrix([1:NumRandomWalkers], 2) + ElementIDinit;
PointMatrix([1:NumRandomWalkers], 3) = PointMatrix([1:NumRandomWalkers], 3) + Lx(ElementIDinit);
PointMatrix([1:NumRandomWalkers], 4) = PointMatrix([1:NumRandomWalkers], 4) + Velocity(ElementIDinit);

Concentration = Element(:, 1) .* 0;
Concentration(ElementIDinit) = 1;

figure(1);
subplot(NumSubPlot, 1, 2)
title('Consentration evolves with time (t = 0)');
xlabel('x (m)'); ylabel('y (m)');
hold on
patch('Vertices', Pnt_y, 'Faces', [1, Element_y(end, [2, 3]), Element_y(1, 4)], 'FaceVertexCData', Pnt_y(:, 1), 'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0); hold on

spo = patch('Vertices', Pnt_y, 'Faces', Element_y(ElementIDinit, :), 'FaceVertexCData', Concentration(ElementIDinit), 'FaceColor', 'flat', 'EdgeAlpha', 0, 'facealpha', 1); hold on
xlim([min(Lx) - 2, max(Lx) + 2])
ylim([-2, 2])
pbaspect([(max(Lx) + 2) * 1, 4, 1]); hold on
caxis([0, 1])
hjk2 = colorbar;
hjk2.Label.String = 'Concentration [-]';

disp('Press enter to continue!')
pause()

VarianceOfX = zeros(NumSteps, 1);
% E_u = zeros(NumSteps, 1);
% E_u_x = zeros(NumSteps, 1);

for i = 1:NumSteps

    AS = find(PointMatrix(:, 1) > 0);

    % ParticleID, ElementID, Position, Velocity
    PointMatrix(AS, 3) = PointMatrix(AS, 3) + deltaT .* PointMatrix(AS, 4) + normrnd(0, 1, [size(AS, 1), 1]) .* sqrt(2 * Dm * deltaT);

    ak1 = find(PointMatrix(AS, 3) >= max(Lx));

    if (isempty(ak1) == 0)
        PointMatrix(AS(ak1), 1) = -1 * PointMatrix(AS(ak1), 1);
        PointMatrix(AS(ak1), 2) = NumPnts - 1;
        PointMatrix(AS(ak1), 3) = max(Lx);
        PointMatrix(AS(ak1), 4) = PointMatrix(AS(ak1), 4) .* 0 + Velocity(NumPnts - 1);
    end

    ak2 = find(PointMatrix(AS, 3) <= min(Lx));

    if (isempty(ak2) == 0)
        PointMatrix(AS(ak2), 1) = -1 * PointMatrix(AS(ak2), 1);
        PointMatrix(AS(ak2), 2) = 1;
        PointMatrix(AS(ak2), 3) = min(Lx);
        PointMatrix(AS(ak2), 4) = PointMatrix(AS(ak2), 4) .* 0 + Velocity(1);
    end

    AF = find(PointMatrix(:, 1) > 0);
    % ParticleID, ElementID, Position, Velocity
    % update elementID
    %     for j = AF'
    %         JKd = find(Lx >= PointMatrix(j, 3), 1);
    %
    %         if (isempty(JKd))
    %             error('error')
    %         else
    %             PointMatrix(j, 2) = JKd - 1;
    %         end
    %
    %     end

    A = repmat([Lx], [1, size(AF, 1)]);
    A = transpose(A);
    
    % ParticleID, ElementID, Position, Velocity
    VarianceOfX(i) = var(PointMatrix(AF, 3));
    % E_u(i) = mean(PointMatrix(AF, 4));
    % E_u_x(i) = mean(PointMatrix(AF, 3) .* PointMatrix(AF, 4));

    cond = bsxfun(@ge, A, PointMatrix(AF, 3));
    [ok, idx] = max(cond, [], 2);

    if (sum(ok) ~= size(ok, 1))
        error('error')
    end

    PointMatrix(AF, 2) = idx - 1;

    CC = PointMatrix(AF, 2);

    edges = unique(CC);
    counts = histc(CC(:), edges);

    Concentration = Concentration .* 0;
    Concentration(edges) = counts;
    Concentration = Concentration ./ NumRandomWalkers;

    %     facealpha_t = Concentration .* 0 + 1;
    %     HK = find(Concentration ~= 0);
    %     facealpha_t(HK) = 1;

    AL = find(Concentration > 0);

    figure(1);
    subplot(NumSubPlot, 1, 2)
    title(['Concentration evolves with time (t = ', num2str(i), ')']);
    delete(spo)
    spo = patch('Vertices', Pnt_y, 'Faces', Element_y(AL, :), 'FaceVertexCData', Concentration(AL), 'FaceColor', 'flat', 'EdgeAlpha', 0, 'facealpha', 1); hold on
    hjk2 = colorbar;
    hjk2.Label.String = 'Concentration [-]';
    caxis([min(Concentration(AL)), max(Concentration(AL))])
    pause(0.1)
end

TimeRange=[1:NumSteps] .* deltaT;
figure(2)
title('Variance vs. Time'); hold on
scatter(TimeRange, VarianceOfX);
hold on
%scatter(TimeRange, 2 .* E_u_x' .*TimeRange + 2 * Dm .* TimeRange - E_u' .^ 2 .* TimeRange);

for i = 1:100
    f = fittype('a*x + b', 'independent', 'x', 'coefficients', {'a', 'b'});
    if (i == 1)
        [cfun, goodness] = fit(TimeRange', VarianceOfX, f, 'startpoint', [(VarianceOfX(end) - VarianceOfX(1))/(TimeRange(end)-TimeRange(1)), 0]);
    else
        [cfun, goodness] = fit(TimeRange', VarianceOfX, f, 'startpoint', [cfun.a, cfun.b]);
    end
    if (goodness.rsquare > 0.98)
        break
    end
end

hold on
plot(TimeRange, TimeRange .* cfun.a + cfun.b, 'r-');

hold on
xlabel('Time'); ylabel('Variance of x')

disp(['-----------------'])
disp('fitting mode is: a*x + b')
disp(['-----------------'])
disp(['goodness.rsquare=', num2str(goodness.rsquare)]);
disp(['a =', num2str(cfun.a),  ', 2 * Dm = ', num2str(2*Dm)]);
%disp(['b =', num2str(cfun.b), ', 2 * Dm = ', num2str(2*Dm)]);
%disp(['c =', num2str(cfun.c)]);
disp(['-----------------'])
