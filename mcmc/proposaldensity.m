function [q,qlog] = proposaldensity(m,p,gamma)
%q = 1./sqrt(2*pi*gamma^2)*exp(-0.5/gamma^2*sum((m(:)-p(:)).^2));
% q = 1;

% parameters:
% m(:,1) thickness 
% m(:,2) vp 
% m(:,3) vs
% m(:,4) density
% m(:,5) anisotropy
% m(:,6) phi
% m(:,7) plunge
% m(:,8) strike
% m(:,9) dip

qlog = 0;
for i = 1:9
    if mean(gamma(:,i)) > 1 % no effective variation
        continue
    end
    if i == 5
        % penalize sum of anisotropy in model
        temp = sum(abs(m(:,i)-p(:,i)));
        gt = mean(gamma(:,i));
        qlog = qlog + -0.5 * sum(((temp).^2) ./ gt.^2 + log(2*pi*gt.^2));
    elseif i == 6 || i == 8
        % angular parameter
        tempa = mod(m(:,i)-p(:,i),1);
        tempb = abs(m(:,i)-p(:,i));
        for j = 1:length(temp)
            temp(j) = min(tempa(j),tempb(j));
        end
        qlog = qlog + -0.5 * sum(((temp).^2) ./ gamma(:,i).^2 + log(2*pi*gamma(:,i).^2));
    else
        qlog = qlog + -0.5 * sum(((m(:,i)-p(:,i)).^2) ./ gamma(:,i).^2 + log(2*pi*gamma(:,i).^2));
    end
end
q = exp(qlog);


end