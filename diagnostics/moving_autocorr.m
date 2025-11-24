function acfSeries = moving_autocorr(x, windowSize, stepSize, lag)
    x = x(:);  % ensure column vector
    N = length(x);
    nWins = floor((N - windowSize) / stepSize) + 1;
    acfSeries = nan(nWins, 2);  % [center iteration, acf value]

    for i = 1:nWins
        startIdx = (i-1) * stepSize + 1;
        endIdx = startIdx + windowSize - 1;
        win = x(startIdx:endIdx);
        win = win - mean(win);
        acf = xcorr(win, lag, 'biased');
        acf = acf / acf(lag+1);  % Normalize by lag 0 value
        acfSeries(i, :) = [startIdx + windowSize/2, acf(lag+1+1)];
    end
end