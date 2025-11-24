function [a]=tritensr(aa,cc,ff,ll,nn,rho)
% Sets up transversely isotropic tensor with horizontal
% symmetry axis.

a = zeros(3,3,3,3);

a(3,3,3,3)=aa./rho;
a(2,2,2,2)=aa./rho;
a(1,1,1,1)=cc./rho;

a(3,3,2,2)=(aa-2.*nn)./rho;
a(2,2,3,3)=(aa-2.*nn)./rho;

a(3,3,1,1)=ff./rho;
a(1,1,3,3)=ff./rho;

a(2,2,1,1)=ff./rho;
a(1,1,2,2)=ff./rho;

a(2,1,2,1)=ll./rho;
a(1,2,1,2)=ll./rho;
a(1,2,2,1)=ll./rho;
a(2,1,1,2)=ll./rho;

a(1,3,1,3)=ll./rho;
a(3,1,3,1)=ll./rho;
a(1,3,3,1)=ll./rho;
a(3,1,1,3)=ll./rho;

a(3,2,3,2)=nn./rho;
a(2,3,2,3)=nn./rho;
a(3,2,2,3)=nn./rho;
a(2,3,3,2)=nn./rho;
end