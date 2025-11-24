function [L,fc,logL] = LeastSquaresLikelyhood_Hkgen_lowA(R,r,Hk,sigR,sigHk,w)
[M,N] = size(R); % 5 × N
% Data misfit term
residual = (R-r).^2 ./ sigR.^2; % element-wise normalized squared residuals
logtermR = log(2 * pi * sigR.^2); % element-wise log of variance term
% Weighted sum over channels (1:5)
rt = zeros(M,1);
for i = 1:M
    rt(i) = w(i) * sum(residual(i,:) + logtermR(i,:))./N;  % sum over time samples (N)
end
% Constraint term (ensure positivity)
temp = max(0, 1 - Hk);
constraint_misfit = w(6) * ((temp^2) / sigHk^2 + log(2*pi*sigHk^2));
% Total negative log-likelihood (normalized)
% disp(num2str(sum(rt)/constraint_misfit))
logL = -0.5 * (sum(rt) + constraint_misfit)/sum(w);%/(M*N+1);
% Misfit per point (optional output)
fc = (sum(rt) + constraint_misfit) / sum(w);% / (M*N+1);
% Likelihood
L = exp(logL);

% rt = sum(((R-r).^2./sigR.^2),2).*w(1:5);
% temp = 1-Hk;
% if temp<0
%     temp = 0;
% end
% tt = ((temp.^2)./sigHk.^2).*w(6);
% % mf2 = sum(h(:).*((a(:)).^2))/sum(h(:));
% % fc = (sum(rt(:))+tt+mf2.*smoothpar.^2)/2/M/sum(w(:));
% fc = (sum(rt(:))./M+tt)/sum(w(:));
% % disp([num2str(mf2.*smoothpar.^2) ' ' num2str(sum(rt(:))+tt)])
% % likelyhood function
% % L = exp(-fc/2)./((((sigR(1)+sigHk)/2).*2.*pi).^(0.5));
% logL = -0.5.*fc - 0.5.*sum(log(2*pi*sigR(:).^2)/M+sigHk.^2);
% L = exp(logL);
end