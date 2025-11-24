function [inc] = get_theo_inc(baz,slow,alpha,beta,rho)
% get inclination for average crustal model
% input parameter
%   baz:    backazimuth for event
%   slow:   ray parameter for event
% optional input:
%   alpha:  Average P-wave velocity for the bulk crust
%   beta:   Average S-wave velocity for the bulk crust
%   rho:    Average density for the bulk crust
% output parameter
%   inc:    inclination of P-wave

% Copyright 2019 F.Link and G.Rümpker

if nargin < 3
    alpha = 6300;
end
if nargin < 4
    beta = 3640;
end
if nargin < 5
    rho = 2850;
end

a = zeros(3,3,3,3);
a(3,3,3,3) = alpha.^2;
a(2,3,2,3) = beta.^2;

p1 = -slow*cosd(baz);
p2 = -slow*sind(baz);

p1r(1)=1;
p1r(2)=0;
p1r(3)=0;
[~,evec] = Hisotroc(a(:,:,:,:),rho,p1,p2);
nu = evec(4:6,4:6);
nd = evec(4:6,1:3);
md = evec(1:3,1:3);
mu = evec(1:3,4:6);
amp = (-mu+md*inv(nd)*nu)*p1r';

P_amp = amp;%(:,1);
inc = atand(sqrt(sum(P_amp(1:2).^2))/abs(P_amp(3)));
end