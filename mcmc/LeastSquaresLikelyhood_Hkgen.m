function [L,fc] = LeastSquaresLikelyhood_Hkgen(R,r,Hk,sigR,sigHk,w)
[M,~] = size(R);
if nargin <6
    w = ones(M+1,1);
end
% cost function
rt = sum(((R-r)./sigR).^2,2).*w(1:5);
tt = (1-Hk)./sigHk.^2.*w(6);
if tt < 0
    tt = 0;
end
fc = (sum(rt(:))+tt)/2/M/sum(w(:));
% likelyhood function
L = exp(-fc/2)./((((sigR(1)+sigHk)/2).*2.*pi).^(0.5));
end