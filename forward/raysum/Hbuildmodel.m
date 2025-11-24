function [a,ar,rot]=buildmodel(rho,alpha_ml,beta,isoflag,pct_a,pct_b,trend,plunge,strike,dip,nlay)
% Set up model for use.

% Parameter eta, combined with percentages of P and S anisotropy,
% is sufficient to obtain all coefficients for a hexagonally-symmetric
% medium. Fixed from Farra et al., 1991

for ilay=1:nlay

    if(isoflag(ilay))
        % Store isotropic coefficients (only 2 are used)
        a(3,3,3,3,ilay) = alpha_ml(ilay).^2.;
        a(2,3,2,3,ilay) = beta(ilay).^2.;
    else
        % Build anisotropic coefficients for a hexagonally-symmetric medium,
        % after Farra et al, 1991
        d_a=alpha_ml(ilay).*pct_a(ilay)./100.;
        d_b=beta(ilay).*pct_b(ilay)./100.;
        aa=rho(ilay).*(alpha_ml(ilay) - d_a./2.).^2.;
        cc=rho(ilay).*(alpha_ml(ilay) + d_a./2.).^2.;
        ll=rho(ilay).*(beta(ilay) + d_b./2.).^2.;
        nn=rho(ilay).*(beta(ilay) - d_b./2.).^2.;
        ac = rho(ilay)*(alpha_ml(ilay)).^2;
%         ff = -ll+sqrt(aa-ll)*sqrt(cc-ll);
        ff = -ll + sqrt((2.*ac).^2 - 2.*ac*(aa + cc + 2.*ll) + (aa + ll)*(cc + ll));
        % disp(num2str(abs(ff-ff2)))
        % eta = ff/(aa - 2.*ll);
        % disp(num2str(eta))
%         eta = 1.03;
%         ff=eta.*(aa-2..*ll);

        % Get tensor with unrotated axes
        [a_temp]=Htritensr(aa,cc,ff,ll,nn,rho(ilay));

        % Rotate axes:
        rot_axis(1,1)=cos(trend(ilay)).*cos(plunge(ilay));
        rot_axis(2,1)=-sin(trend(ilay));
        rot_axis(3,1)=-cos(trend(ilay)).*sin(plunge(ilay));
        rot_axis(1,2)=sin(trend(ilay)).*cos(plunge(ilay));
        rot_axis(2,2)=cos(trend(ilay));
        rot_axis(3,2)=-sin(trend(ilay)).*sin(plunge(ilay));
        rot_axis(1,3)=sin(plunge(ilay));
        rot_axis(2,3)=0.;
        rot_axis(3,3)=cos(plunge(ilay));

        a(:,:,:,:,ilay)=Hrot_tensor(a_temp,rot_axis);

    end
end


% Make rotator-matrix list:
rot=Hmake_rotator(strike,dip,nlay);

ar = zeros(3,3,3,3,nlay,2);
% Pre-rotate tensors (skip isotropic ones)
for ilay=1:nlay
    if(~(isoflag(ilay)))
        %            Upper interface:
        ar(:,:,:,:,ilay,2)=Hrot_tensor(a(:,:,:,:,ilay), rot(:,:,ilay));
        %            Lower interface. Bottom layer is a half-space.
        if(ilay < nlay)
            ar(:,:,:,:,ilay,1)=Hrot_tensor(a(:,:,:,:,ilay), rot(:,:,ilay+1));
        end
    end
end


end