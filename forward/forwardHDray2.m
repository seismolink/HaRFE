function [hdc,rt,tt,wr,wt,vt] = forwardHDray2(model,param,geom,t_ref,phaselist,nseg)

M = size(phaselist);
numph = M(3);

Nsample2 = length(t_ref{1});

myData = param;
Nsample = myData(2);
dt = myData(3);
width = myData(4);
identrot = myData(7);
maxt = myData(10);

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
vpref = (max(vp)+200)/1000;
pray(pray*1000>1/vpref) = 1/vpref/1000;

% Define gaussian filter
df = 1/(Nsample*dt);
do=2*pi*df;
oo=(0:Nsample-1)*do;
a2=4*(width*2*pi).^2;
centre_gauss = 0;
gw = exp(-((oo-centre_gauss*pi*2).*(oo-centre_gauss*pi*2))/a2);

bootg = zeros(10,10,Nsample2);
bootd = zeros(10,Nsample2);
Wn = bootd;
Wd = bootd;
tt = zeros(length(baz),Nsample);
rt = zeros(length(baz),Nsample);
vt = zeros(length(baz),Nsample);
wr = rt;
wt = tt;
rt0 = rt;
vt0 = rt;
tt0 = tt;

time = linspace(0,(param(2)-1).*param(3),param(2))+param(6);
index = 1:length(time);
uflag = zeros(size(baz));

pp = gcp('nocreate');
if isempty(pp)
    if isempty(getenv('SLURM_CPUS_ON_NODE'))
        nWorkers = feature('NumCores');
    else
        nWorkers = str2double(getenv('SLURM_CPUS_ON_NODE'));
    end
    pp = parpool(nWorkers);
end
parfor ii = 1:length(baz)
% if isscalar(baz)
%     hdc = 0;
%     continue
% end
% for ii = 1:length(baz)

% Set up model for calculation, make rotators
% -------------------------------------------
% Variables
bazmod = baz(ii)/180*pi;
sta_dx = 0;
sta_dy = 0;


% Perform calculation for amplitudes and arrival times
% ----------------------------------------------------
amp_in = 1;
[travel_time, amplitude] = Hget_arrivals(H,rho,avp==0,strike,dip,aa,ar_list,rot,bazmod,pray(ii),sta_dx,sta_dy,phaselist,1,nseg,numph,nlay,amp_in);
[amplitude]=Hnorm_arrivals(amplitude,bazmod,pray(ii),vp(1),vs(1),rho(1),1,numph, 1,1);

V = -amplitude(3,:); % Z
[~,I] = max(abs(V));
H1 = amplitude(1,:)./V(I); % NS
H2 = amplitude(2,:)./V(I); % EW
V = V./V(I);

V(travel_time==0) = [];
H1(travel_time==0) = [];
H2(travel_time==0) = [];
travel_time(travel_time==0) = [];
if length(travel_time) < 1
    disp('Error in forward calculation! No travel time calculated.')
    continue
end
traveltime = travel_time-travel_time(1);%+shift;

% Prepare traces
% --------------
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
    % disp(num2str(inc))
    inc = -inc;
end
trVa = trV .* cosd(inc) + trH1 .* sind(inc) .* cosd(baz(ii)) + trH2 .* sind(inc) .* sind(baz(ii));
trH1a = trV .* sind(inc) - trH1 .* cosd(inc) .* cosd(baz(ii)) - trH2 .* cosd(inc) .* sind(baz(ii));
trH2a = -(trH1 .* sind(baz(ii)) - trH2 .* cosd(baz(ii)));

[~,I] = max(abs(trVa));
fnorm = trVa(I);
V = fft(trVa./fnorm);
H1 = fft(trH1a./fnorm);
H2 = fft(trH2a./fnorm);

temp = real(ifft(H1.*(gw)));
rt(ii,:) = temp;
% rt(ii,1:length(temp)) = temp;
temp = real(ifft(H2.*(gw)));
tt(ii,:) = temp;
% tt(ii,1:length(temp)) = temp;
temp = real(ifft(V.*gw));
vt(ii,:) = temp;
% vt(ii,1:length(temp)) = temp;

rt0(ii,:) = rt(ii,:);
vt0(ii,:) = vt(ii,:);
tt0(ii,:) = tt(ii,:);
vt(ii,:) = vt(ii,:)./max(abs(rt(ii,:)));
tt(ii,:) = tt(ii,:)./max(abs(rt(ii,:)));
rt(ii,:) = rt(ii,:)./max(abs(rt(ii,:)));

rat = 0.1;
[~,idx] = findpeaks(abs(rt(ii,:)));
temp = rt(ii,:);
ar = temp(idx);
temp = wr(ii,:);
temp(idx(abs(ar)>rat.*max(abs(ar)))) = sign(ar(abs(ar)>rat.*max(abs(ar))));
wr(ii,:) = temp;
% wr(ii,idx(abs(ar)>rat.*max(abs(ar)))) = sign(ar(abs(ar)>rat.*max(abs(ar))));
[~,idxt] = findpeaks(abs(tt(ii,:)));
temp = tt(ii,:);
at = temp(idxt);
temp = wt(ii,:);
temp(idxt(abs(at)>rat.*max(abs(ar)))) = sign(at(abs(at)>rat.*max(abs(ar))));
wt(ii,:) = temp;
% wt(ii,idxt(abs(at)>rat.*max(abs(ar)))) = sign(at(abs(at)>rat.*max(abs(ar))));
wr(ii,:) = real(ifft(fft(wr(ii,:)).*gw));
wt(ii,:) = real(ifft(fft(wt(ii,:)).*gw));
temp = wr(ii,:);
fnormr = max(abs(wr(ii,:)-temp(1)));
if fnormr == 0
    fnormr = 1;
end
wr(ii,:) = (wr(ii,:)-temp(1))./fnormr;
temp = wt(ii,:);
fnormt = max(abs(wt(ii,:)-temp(1)));
if fnormt == 0
    fnormt = 1;
end
wt(ii,:) = (wt(ii,:)-temp(1))./fnormt;

uflag(ii) = 1;
end


for ii = 1:length(baz)
    if uflag(ii)==0
        continue
    end
    if isscalar(baz)
        hdc = 0;
        continue
    end
% ATTENTION!!! Harmonic decomposition requires right handed coordinate
% system ->> change the sign on the transverse component T=-T

H1 = rt(ii,:)'./max(abs(rt(ii,:)));
H2 = -tt(ii,:)'./max(abs(rt(ii,:)));
    
    bazb = baz(ii);
    cbaz = cosd(bazb);
    sbaz = sind(bazb);
    c2baz = cbaz*cbaz-sbaz*sbaz;
    s2baz = 2.0*cbaz*sbaz;
    
    bf(1) = 1.0;
    % bf(1) = 3.0;
    bf(2) = cbaz;
    bf(3) = sbaz;
    bf(4) = c2baz;
    bf(5) = s2baz;
    bf(6) = 0.0;
    bf(7) = cbaz;
    bf(8) = sbaz;
    bf(9) = c2baz;
    bf(10) = s2baz;
    
    bft(1) = 0.0;
    % bft(1) = 0.3;
    bft(2) = sbaz;
    bft(3) = -cbaz;
    bft(4) = s2baz;
    bft(5) = -c2baz;
    bft(6) = 1.0;
    bft(7) = -sbaz;
    bft(8) = cbaz;
    bft(9) = -s2baz;
    bft(10) = c2baz;
    
    
    bfdyad = bf'*bf;
    bftdyad = bft'*bft;
    
    rff = interp1(time,H1,t_ref{ii})';
    rfft = -interp1(time,H2,t_ref{ii})';
    rff(isnan(rff)) = 0;
    rfft(isnan(rfft)) = 0;
    % if ii > 9
    %     rff = rff*10;
    %     rfft = rfft*10;
    % end

    bootg = bootg + repmat(bfdyad,1,1,Nsample2) + repmat(bftdyad,1,1,Nsample2);
    bootd = bootd + bf'*rff' + bft'*rfft';

    % Wn = Wn + bf'*rff' + bft'*rfft';
    % Wd = Wd + abs(bf'*rff' + bft'*rfft');
end
rt = rt0;
tt = tt0;
vt = vt0;
rt(:,time>maxt) = [];
tt(:,time>maxt) = [];
vt(:,time>maxt) = [];
wt(:,time>maxt) = [];
wr(:,time>maxt) = [];
% bootd = bootd.*(abs(Wn)./Wd);
for i = 1:Nsample2
    bootm(i,:) = bootg(:,:,i)\bootd(:,i);
end
% for i = 1:10
%     hdc(:,i) = interp1(t_ref{end},bootm(:,i),time);%(real(ifft(bootm(:,i).*gw')));
%     % hdc(:,i) = (real(ifft(bootm(:,i).*gw')));
% end
% hdc(isnan(hdc)) = 0;
% temp = abs(hdc(time<=param(8),:));
fn = max(abs(bootm(:)));
hdc = bootm./fn;
% hdc(time>param(8),:) = 0;
% hdc(time>param(8),:) = [];
% time(time>param(8)) = [];
% hdc(time<param(6)) = [];
end