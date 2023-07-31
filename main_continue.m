clc
close all

% run more steps
NumSteps_continue = 200;

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
figure(1);
subplot(NumSubPlot, 1, 2)
title(['Concentration evolves with time (step = ', num2str(NumSteps), ')']);
xlabel('x (m)'); ylabel('y (m)');
hold on
patch('Vertices', Pnt_y, 'Faces', [1, Element_y(end, [2, 3]), Element_y(1, 4)], 'FaceVertexCData', Pnt_y(:, 1), 'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0); hold on

% spo = patch('Vertices', Pnt_y, 'Faces', Element_y(ElementIDinit, :), 'FaceVertexCData', Concentration(ElementIDinit), 'FaceColor', 'flat', 'EdgeAlpha', 0, 'facealpha', 1); hold on
xlim([min(Lx) - 2, max(Lx) + 2])
ylim([-2, 2])
pbaspect([(max(Lx) + 2) * 1, 4, 1]); hold on
caxis([0, 1])
hjk2 = colorbar;
hjk2.Label.String = 'Concentration [-]';

% enlarge the vectors
VarianceOfX = [VarianceOfX; zeros(NumSteps_continue, 1)];
E_u = [E_u; zeros(NumSteps_continue, 1)];
E_u_x = [E_u_x; zeros(NumSteps_continue, 1)];

for i = NumSteps + 1:NumSteps + NumSteps_continue

    AS = find(RandomWalkerMatrix(:, 1) > 0);

    % ParticleID, ElementID, Position, Velocity
    RandomWalkerMatrix(AS, 3) = RandomWalkerMatrix(AS, 3) + deltaT .* RandomWalkerMatrix(AS, 4) + normrnd(0, 1, [size(AS, 1), 1]) .* sqrt(2 * Dm * deltaT);

    ak1 = find(RandomWalkerMatrix(AS, 3) >= max(Lx));

    if (isempty(ak1) == 0)
        RandomWalkerMatrix(AS(ak1), 1) = -1 * RandomWalkerMatrix(AS(ak1), 1);
        RandomWalkerMatrix(AS(ak1), 2) = NumPnts - 1;
        RandomWalkerMatrix(AS(ak1), 3) = max(Lx);
        RandomWalkerMatrix(AS(ak1), 4) = RandomWalkerMatrix(AS(ak1), 4) .* 0 + Velocity(NumPnts - 1);
    end

    ak2 = find(RandomWalkerMatrix(AS, 3) <= min(Lx));

    if (isempty(ak2) == 0)
        RandomWalkerMatrix(AS(ak2), 1) = -1 * RandomWalkerMatrix(AS(ak2), 1);
        RandomWalkerMatrix(AS(ak2), 2) = 1;
        RandomWalkerMatrix(AS(ak2), 3) = min(Lx);
        RandomWalkerMatrix(AS(ak2), 4) = RandomWalkerMatrix(AS(ak2), 4) .* 0 + Velocity(1);
    end

    AF = find(RandomWalkerMatrix(:, 1) > 0);
    A = repmat([Lx], [1, size(AF, 1)]);
    A = transpose(A);
    % ParticleID, ElementID, Position, Velocity
    cond = bsxfun(@ge, A, RandomWalkerMatrix(AF, 3));
    [ok, idx] = max(cond, [], 2);
    clear A

    if (sum(ok) ~= size(ok, 1))
        error('error')
    end

    % update ElementID and velocity
    RandomWalkerMatrix(AF, 2) = idx - 1;
    RandomWalkerMatrix(AF, 4) = Velocity(RandomWalkerMatrix(AF, 2));

    % record expectation of u and x
    VarianceOfX(i) = var(RandomWalkerMatrix(AF, 3));
    E_u(i) = mean(RandomWalkerMatrix(AF, 4));
    E_u_x(i) = mean(RandomWalkerMatrix(AF, 3) .* RandomWalkerMatrix(AF, 4));
    % E_x(i) = mean(RandomWalkerMatrix(AF, 3));

    % ElementID
    CC = RandomWalkerMatrix(AF, 2);

    edges = unique(CC);
    counts = histc(CC(:), edges);

    Concentration = Concentration .* 0;
    Concentration(edges) = counts;
    Concentration = Concentration ./ NumRandomWalkers;

    AL = find(Concentration > 0);

    figure(1);
    subplot(NumSubPlot, 1, 2)
    title(['Concentration evolves with time (step = ', num2str(i), ')']);
    delete(spo)
    spo = patch('Vertices', Pnt_y, 'Faces', Element_y(AL, :), 'FaceVertexCData', Concentration(AL), 'FaceColor', 'flat', 'EdgeAlpha', 0, 'facealpha', 1); hold on
    hjk2 = colorbar;
    hjk2.Label.String = 'Concentration [-]';
    caxis([min(Concentration(AL)), max(Concentration(AL))])
    pause(0.1)
end

NumSteps = NumSteps + NumSteps_continue;

TimeRange = [1:NumSteps] .* deltaT;
figure(2)
title('Variance vs. Time'); hold on
xlabel('Time'); ylabel('Variance of x'); hold on
HYU1 = scatter(TimeRange, VarianceOfX);
hold on
%scatter(TimeRange, 2 .* E_u_x' .*TimeRange + 2 * Dm .* TimeRange - E_u' .^ 2 .* TimeRange);

if (IfLinearFitting == 1)
    InitTimePnt = 10;
    TimeRange_II = TimeRange([InitTimePnt:end]);
    VarianceOfX_II = VarianceOfX([InitTimePnt:end]);

    for i = 1:100
        f = fittype('a*x + b', 'independent', 'x', 'coefficients', {'a', 'b'});

        if (i == 1)
            [cfun, goodness] = fit(TimeRange_II', VarianceOfX_II, f, 'startpoint', [(VarianceOfX_II(end) - VarianceOfX_II(1)) / (TimeRange_II(end) - TimeRange_II(1)), 0]);
        else
            [cfun, goodness] = fit(TimeRange_II', VarianceOfX_II, f, 'startpoint', [cfun.a, cfun.b]);
        end

        if (goodness.rsquare > 0.98)
            break
        end

    end

    hold on
    HYU2 = plot(TimeRange_II, TimeRange_II .* cfun.a + cfun.b, 'r-', 'linewidth', 2);
    hold on
    xlabel('Time'); ylabel('Variance of x')
    disp(['-----------------'])
    disp('fitting mode is: a*x + b')
    disp(['-----------------'])
    disp(['goodness.rsquare=', num2str(goodness.rsquare)]);
    disp(['a =', num2str(cfun.a), ', 2 * Dm = ', num2str(2 * Dm)]);
    disp(['b =', num2str(cfun.b)]);
    %disp(['c =', num2str(cfun.c)]);
    disp(['-----------------'])
    legend([HYU1 HYU2], {'Variance vs. time', 'model: 2 Dm t'});
else
    % plot(TimeRange, 2 .* TimeRange' .* E_u_x + 2 * Dm .* TimeRange'  - E_u .^ 2 .* TimeRange' .^ 2, 'r-', 'linewidth', 2);
    %E_u_x_Integrated = zeros(NumSteps, 1);
    %E_u_Integrated = zeros(NumSteps, 1);
    %E_x_Integrated = E_x;

    E_u_x_Integrated = cumtrapz([1:NumSteps] .* deltaT, E_u_x);
    E_u_Integrated = cumtrapz([1:NumSteps] .* deltaT, E_u);

    % HYU2 = plot(TimeRange, TimeRange' .* E_u_x + 2 * Dm .* TimeRange' - E_u .^ 2 .* TimeRange' .^ 2, 'r-', 'linewidth', 2);
    HYU2 = plot(TimeRange, 2 .* E_u_x_Integrated + 2 * Dm .* TimeRange' - E_u_Integrated .^ 2, 'r-', 'linewidth', 2);
    % HYU2 = plot(TimeRange, TimeRange' .* E_x_Integrated .* mean(Velocity) + 2 * Dm .* TimeRange'  -  mean(Velocity) .^ 2 .* TimeRange' .^ 2, 'r-', 'linewidth', 2);
    legend([HYU1 HYU2], {'Variance vs. time', 'model: the completed form'})
end
