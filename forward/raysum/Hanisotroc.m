function [eval,evec]=anisotroc(a,rho,p1,p2)
%   Obtain anisotropic eigenvectors/values. Rather hairy.

%       Build partion matrices CIJ where CIJ(k,l)=rho*a(k,i,l,j)
cc = rho.*permute(a,[1 3 2 4]);

%       T=(-p1*C13-p2*C23)*iC33
t = (-p1*cc(:,:,1,3)-p2*cc(:,:,2,3))/cc(:,:,3,3);

%       S=rho*eye(3,3)-sum(i=1:2,j=1:2,p(i)*p(j)*(Cij-Ci3*C33^(-1)*C3j))
p(1)=p1;
p(2)=p2;
s = rho.*eye(3)-(p1*p1*(cc(:,:,1,1)-cc(:,:,1,3)*(cc(:,:,3,3)\cc(:,:,3,1)))+...
    p1*p2*(cc(:,:,1,2)-cc(:,:,1,3)*(cc(:,:,3,3)\cc(:,:,3,2)))+...
    p2*p1*(cc(:,:,2,1)-cc(:,:,2,3)*(cc(:,:,3,3)\cc(:,:,3,1)))+...
    p2*p2*(cc(:,:,2,2)-cc(:,:,2,3)*(cc(:,:,3,3)\cc(:,:,3,2))));

%       Build system matrix AA=[T',iC33;S,T]
aa(1:3,1:3) = t';
aa(1:3,4:6) = inv(cc(:,:,3,3));
aa(4:6,1:3) = s;
aa(4:6,4:6) = t;

%       Obtain eigenvalues/vectors using EISPACK.
[V,D] = eig(aa);
d = diag(D);
evalr = real(d);
evali = imag(d);
evecr = real(V);
eveci = imag(V);

%        Sort evecs/evals
[eval,evec]=Hsort_evec(evalr,evali,evecr,eveci);

%        Normalize evecs to unit displacement
for j=1:6
    xnorm=sqrt(real(evec(1,j)).^2+real(evec(2,j)).^2+ real(evec(3,j)).^2);
    evec(:,j)=evec(:,j)./xnorm;
end

end