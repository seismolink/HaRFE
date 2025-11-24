function [eval,evec]=isotroc(a,rho,p1,p2)
%   Obtain isotropic eigenvectors/values. This is a little obscure,
%   but it's just the Zoeppritz equations.

%   Load up some convenient values
vp2=a(3,3,3,3);
vs2=a(2,3,2,3);
mu=rho.*vs2;
pp=p1.*p1+p2.*p2;

%    Set eigenvalues
qdp=sqrt(1./vp2-pp);
qds=sqrt(1./vs2-pp);
qup=-qdp;
qus=-qds;
eval(1)=qdp;
eval(2)=qds;
eval(3)=qds;
eval(4)=qup;
eval(5)=qus;
eval(6)=qus;

%    Set eigenvector matrix, column by column,
%      N(1:6,1)=[p1,p2,qdp,2*mu*p1*qdp,2*mu*p2*qdp,(rho-2*mu*pp)]'
evec(1,1)=p1;
evec(2,1)=p2;
evec(3,1)=qdp;
evec(4,1)=2..*mu.*p1.*qdp;
evec(5,1)=2..*mu.*p2.*qdp;
evec(6,1)=rho-2..*mu.*pp;
%      N(1:6,2)=[p1,p2,-pp/qds,p1*N(6,1)/qds,p2*N(6,1)/qds,-2*mu*pp]';
evec(1,2)=p1;
evec(2,2)=p2;
evec(3,2)=-pp./qds;
evec(4,2)=p1.*evec(6,1)./qds;
evec(5,2)=p2.*evec(6,1)./qds;
evec(6,2)=-2..*mu.*pp;
%      N(1:6,3)=[-p2,p1,0,-p2*qds*mu,p1*qds*mu,0]'
evec(1,3)=-p2;
evec(2,3)=p1;
evec(3,3)=0.;
evec(4,3)=-p2.*qds.*mu;
evec(5,3)=p1.*qds.*mu;
evec(6,3)=0.;
%      N(1:6,4)=[p1,p2,qup,2*mu*p1*qup,2*mu*p2*qup,N(6,1)]';
evec(1,4)=p1;
evec(2,4)=p2;
evec(3,4)=qup;
evec(4,4)=2..*mu.*p1.*qup;
evec(5,4)=2..*mu.*p2.*qup;
evec(6,4)=evec(6,1);
%      N(1:6,5)=[p1,p2,-pp/qus,p1*N(6,1)/qus,p2*N(6,1)/qus,-2*mu*pp]';
evec(1,5)=p1;
evec(2,5)=p2;
evec(3,5)=-pp./qus;
evec(4,5)=p1.*evec(6,1)./qus;
evec(5,5)=p2.*evec(6,1)./qus;
evec(6,5)=-2..*mu.*pp;
%      N(1:6,6)=[-p2,p1,0,-p2*qus*mu,p1*qus*mu,0]'
evec(1,6)=-p2;
evec(2,6)=p1;
evec(3,6)=0.;
evec(4,6)=-p2.*qus.*mu;
evec(5,6)=p1.*qus.*mu;
evec(6,6)=0.;

%      Normalize wrt displacement magnitude:
for j=1:6
    xnorm=sqrt(real(evec(1,j)).^2+real(evec(2,j)).^2+ real(evec(3,j)).^2);
    evec(:,j) = evec(:,j)./xnorm;
end

end