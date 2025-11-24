function [acfSeries, centers] = moving_autocorr_multilag(x, windowSize, stepSize, maxLag)
    % Inputs:
    %   x          - MCMC chain (vector)
    %   windowSize - number of steps per window
    %   stepSize   - how far to slide the window
    %   maxLag     - number of autocorrelation lags to compute

    x = x(:);  % ensure column
    N = length(x);
    nWins = floor((N - windowSize) / stepSize) + 1;
    acfSeries = nan(nWins, maxLag+1);  % rows = windows, cols = lags
    centers = nan(nWins, 1);

    for i = 1:nWins
        startIdx = (i-1) * stepSize + 1;
        endIdx = startIdx + windowSize - 1;
        win = x(startIdx:endIdx);
        win = win - mean(win);
        [acf_full, lags] = xcorr(win, maxLag, 'biased');
        acf = acf_full(lags >= 0);
        acf = acf / acf(1);  % normalize by lag-0 value

        acfSeries(i, :) = acf(:)';
        centers(i) = startIdx + windowSize / 2;
    end
end