clear all
close all

% nSteps = 10000;
% nDims = 2;
% rho = 0.98;
% chain = zeros(nSteps, nDims);
% chain(1,:) = randn(1, nDims);
% for t = 2:nSteps
%     chain(t,:) = rho * chain(t-1,:) + sqrt(1 - rho^2) * randn(1, nDims);
% end
addname = 'Panama3';
statname = 'BCIP';
chainname = 'a_5';

load(['D:\Panama\Inversions\Panama3\BCIP\MCMCresult' chainname '.mat'])

% Parameters
maxLag = 1000;  % for full and moving ACF
windowSize = 100;
stepSize = 10;
lagsToPlot = [1, 3, 5, 15, 30, 50];

nn = size(MC.post.para);
chain = reshape(MC.post.para,nn(1)*nn(2),nn(3))';
nDims = nn(1)*nn(2);

% chain = MC.post.full.stat(1:end,:);
% nn = size(chain);
% nDims = nn(2);

% Diagnostics
for d = 1:1%nDims
    x = chain(:,d);
    x_demeaned = x - mean(x);

    % --- Subplot block for parameter d ---
    fig = figure('Name', ['Diagnostics for Param ', num2str(d)], 'Position', [100,100,1400,500]);

    % (1) Trace
    subplot(1,3,1);
    plot(x);
    xlabel('Iteration'); ylabel('Value');
    title(['Trace: Param ', num2str(d)]);

    % (2) Full-chain ACF
    [acf_full, lags] = xcorr(x_demeaned, maxLag, 'biased');
    acf = acf_full(lags >= 0);
    acf = acf / acf(1);
    subplot(1,3,2);
    stem(0:maxLag, acf, 'filled');
    xlabel('Lag'); ylabel('Autocorr');
    title(['ACF (Full chain) - Param ', num2str(d)]);
    ylim([-0.1 1]);

    % (3) Moving ACF at multiple lags
    [acfTrack, centers] = moving_autocorr_multilag(x, windowSize, stepSize, maxLag);
    subplot(1,3,3);
    hold on;
    colors = lines(length(lagsToPlot));
    for i = 1:length(lagsToPlot)
        lagIdx = lagsToPlot(i) + 1;
        plot(centers, acfTrack(:, lagIdx), '-o', 'Color', colors(i,:), ...
             'DisplayName', ['Lag ', num2str(lagsToPlot(i))]);
    end
    xlabel('Iteration (center of window)');
    ylabel('Autocorrelation');
    title(['Moving ACF - Param ', num2str(d)]);
    legend show;
    grid on;
end
print(fig,['Diagnostics' addname '_' statname '_' chainname '_param' num2str(d) '.png'],'-r300','-dpng')