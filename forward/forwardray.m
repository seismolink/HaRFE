function [rt,tt,wr,wt] = forwardray(model,param,geom,phaselist,nseg)

M = size(phaselist);
numph = M(3);

myData = param;
Nsample = myData(2);
dt = myData(3);
width = myData(4);
identrot = myData(7);

myData = model;
H = myData(:,1)';
rho = myData(:,2)';
vp = myData(:,3)';
vs = myData(:,4)';
avp = myData(:,6)';
avs = myData(:,7)';
phi = myData(:,8)'/180*pi;
plunge = myData(:,9)'/180*pi;
strike = myData(:,10)'/180*pi;
dip = myData(:,11)'/180*pi;
nlay = length(H);
[aa,ar_list,rot] = Hbuildmodel(rho,vp,vs,avp==0,avp,avs,phi,plunge,strike,dip,nlay);

myData = geom;
baz = myData(:,1)';
pray = myData(:,2)./1000';

% Define gaussian filter
df = 1/(Nsample*dt);
do=2*pi*df;
oo=(0:Nsample-1)*do;
a2=4*(width*2*pi).^2;
centre_gauss = 0;
gw = exp(-((oo-centre_gauss*pi*2).*(oo-centre_gauss*pi*2))/a2);

tt = zeros(length(baz),Nsample);
rt = zeros(length(baz),Nsample);
wr = rt;
wt = tt;
for ii = 1:length(baz)
    
% Set up model for calculation, make rotators
% -------------------------------------------
% Variables
bazmod = baz(ii)/180*pi;
sta_dx = 0;
sta_dy = 0;


% Perform calculation for amplitudes and arrival times
% ----------------------------------------------------
amp_in = 1;
% if pray(ii)>0.1
%     pr = 0.1;
% else
    pr = pray(ii);
    % bazmod = round(baz(ii))/180*pi;
% end
[travel_time, amplitude] = Hget_arrivals(H,rho,avp==0,strike,dip,aa,ar_list,rot,bazmod,pr,sta_dx,sta_dy,phaselist,1,nseg,numph,nlay,amp_in);
[amplitude]=Hnorm_arrivals(amplitude,bazmod,pr,vp(1),vs(1),rho(1),1,numph, 1,1);

V = -amplitude(3,:); % Z
[~,I] = max(abs(V));
H1 = amplitude(1,:)./V(I); % NS
H2 = amplitude(2,:)./V(I); % EW
V = V./V(I);

V(travel_time==0) = [];
H1(travel_time==0) = [];
H2(travel_time==0) = [];
travel_time(travel_time==0) = [];
traveltime = travel_time-travel_time(1);%+shift;

% Prepare traces
% --------------
time = linspace(0,(param(2)-1).*param(3),param(2))+param(6);
index = 1:length(time);
postime = round(interp1(time,index,traveltime));

indx = unique(postime);
ampV = zeros(size(indx));
ampH1 = zeros(size(indx));
ampH2 = zeros(size(indx));
for j = 1:length(indx)
    ampV(j) = sum(V(indx(j)==postime));
    ampH1(j) = sum(H1(indx(j)==postime));
    ampH2(j) = sum(H2(indx(j)==postime));
end
trV = zeros(size(time));
trV(indx) = ampV;
trH1 = zeros(size(time));
trH1(indx) = ampH1;
trH2 = zeros(size(time));
trH2(indx) = ampH2;

if identrot == 1
    inc = 0;
else
    [inc] = get_theo_inc(baz(ii),pray(ii),vp(1),vs(1),rho(1));
    disp(num2str(inc))
    inc = -inc;
end
trVa = trV .* cosd(inc) + trH1 .* sind(inc) .* cosd(baz(ii)) + trH2 .* sind(inc) .* sind(baz(ii));
trH1a = trV .* sind(inc) - trH1 .* cosd(inc) .* cosd(baz(ii)) - trH2 .* cosd(inc) .* sind(baz(ii));
trH2a = (trH1 .* sind(baz(ii)) - trH2 .* cosd(baz(ii)));
% trVa = trV .* cosd(inc) - trH1 .* sind(inc) .* cosd(baz(ii)) - trH2 .* sind(inc) .* sind(baz(ii));
% trH1a = trV .* sind(inc) + trH1 .* cosd(inc) .* cosd(baz(ii)) + trH2 .* cosd(inc) .* sind(baz(ii));
% trH2a = -(trH1 .* sind(baz(ii)) + trH2 .* cosd(baz(ii)));

[~,I] = max(abs(trVa));
fnorm = trVa(I);
V = fft(trVa./fnorm);
H1 = fft(trH1a./fnorm);
H2 = fft(trH2a./fnorm);

temp = real(ifft(H1.*(gw)));
rt(ii,1:length(temp)) = temp;
temp = real(ifft(H2.*(gw)));
tt(ii,1:length(temp)) = temp;
% temp = real(ifft(V.*gw));
% vt(ii,1:length(temp)) = temp;

% vt(ii,:) = vt(ii,:)./max(abs(rt(ii,:)));
tt(ii,:) = tt(ii,:)./max(abs(rt(ii,:)));
rt(ii,:) = rt(ii,:)./max(abs(rt(ii,:)));

rat = 0.1;
[~,idx] = findpeaks(abs(rt(ii,:)));
ar = rt(ii,idx);
wr(ii,idx(abs(ar)>rat.*max(abs(ar)))) = sign(ar(abs(ar)>rat.*max(abs(ar))));
[~,idxt] = findpeaks(abs(tt(ii,:)));
at = tt(ii,idxt);
wt(ii,idxt(abs(at)>rat.*max(abs(ar)))) = sign(at(abs(at)>rat.*max(abs(ar))));
wr(ii,:) = real(ifft(fft(wr(ii,:)).*gw));
wt(ii,:) = real(ifft(fft(wt(ii,:)).*gw));
wr(ii,:) = (wr(ii,:)-wr(ii,1))./max(abs(wr(ii,:)-wr(ii,1)));
wt(ii,:) = (wt(ii,:)-wt(ii,1))./max(abs(wt(ii,:)-wt(ii,1)));

end
