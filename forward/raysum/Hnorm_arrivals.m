function [amp]=norm_arrivals(amp,baz,slow,alpha_ml,beta,rho,ntr,numph, arr,comp)
%      Normalize arrival amplitudes amp by the amplitude of arrival arr,
%      component (in P-SV-SH system) comp. The arrivals are assumed to
%      be in NS-EW-Z coordinates.

a(3,3,3,3)=alpha_ml.^2;
a(2,3,2,3)=beta.^2;

for itr=1:ntr

p1=-slow(itr).*cos(baz(itr));
p2=-slow(itr).*sin(baz(itr));
[~,evec] = Hisotroc(a,rho,p1,p2);
md = evec(1:3,1:3);
mu = evec(1:3,4:6);
nd = evec(4:6,1:3);
nu = evec(4:6,4:6);
% invop = inv(complex(1,0)*mu+complex(-1,0)*nu/nd*md);
% 
% temp1 = inv(nd);
% temp2 = nu*temp1;
% temp3 = md*temp2;
% temp4 = complex(1,0).*mu+complex(-1,0).*temp3;
% invop = inv(temp4);

% invop = inv(complex(1,0).*mu+complex(-1,0).*(md*nu/nd));

% [nd,invop]=cmatinv3(nd,invop);
% [invop,nu,wrk]=cmatmul3(invop,nu,wrk);
% [md,wrk,invop]=cmatmul3(md,wrk,invop);
% [dumvar1,mu,dumvar3,invop,wrk]=cmatlincomb3(complex(1.,0),mu,complex(-1.,0),invop,wrk);
% [wrk,invop]=cmatinv3(wrk,invop);

% invop = inv(mu-md*inv(nd)*nu);


% for j=1:3
%     u(j)=amp(j,arr,itr);
% end
u = squeeze(amp(:,arr,itr));
% keyboard
% v = u'*invop;
v =  u'/(mu-md*(nd\nu))';
normamp=-real(v(comp));

% disp(num2str(u(comp)));
% disp(num2str(-real(v(comp))));

amp(:,:,itr)=amp(:,:,itr)./normamp;
% for iph=1:numph
% for j=1:3
% % if(normamp > 0.)
% amp(j,iph,itr)=amp(j,iph,itr)./normamp;
% % else
% % amp(j,iph,itr)=0.;
% % end
% end
% end

end


end