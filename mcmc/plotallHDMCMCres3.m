function plotallHDMCMCres3(sel_data)
close all
if nargin < 1
    inputflag = 0;
else
    inputflag = 1;
end
% Bayesian inversion for anisotropic crustal model allowing for dipping
% layers based on a Markov chain Monte Carlo algorithm with Metropolis
% Hastings acceptance criterion
% ---------------------------
% Data adaptive version

NN1 = 5000;
NN2 = 200;

% get current folder name
cur_dir = pwd;
cd ..
programpath = pwd;
[~, df, ~] = fileparts(programpath); 

% add all sub directories to search path
addpath(genpath(['../',df]));
cd(cur_dir)

if ~inputflag
    results_dir = input('Define path to results: ','s');
    input_file = input('Define absolute path to input file: ','s');
    modpath = input('Path to modelfiles as "stationname.txt" for comparison: ','s');
    geompath = input('Path to geometry file: ','s');
    bfflag = input('Use mean (0) or best fit (1) as result? : ');
else
    results_dir = sel_data.results_dir;
    input_file = sel_data.input_file;
    modpath = sel_data.modpath;
    bfflag = sel_data.bfflag;
    geompath = '';
end

savepath = results_dir;
if ~isempty(modpath)
    addname = 'Final_mod';
else
    addname = 'Final';
end
if bfflag
    addname = [addname 'bf'];
end

A = dir([results_dir '/MCMCresult*.mat']);
if isempty(A)
    return
end
n = 1;
for i = 1:length(A)
    load([A(i).folder '/' A(i).name])
    for j = 1:length(MC)
        MCa(n) = MC(j);
        n = n+1;
    end
end
MC = MCa;
Nmc = length(MC);
N = MC(1).N;
param = MC(1).param;
w = MC(1).w;
b = MC(1).b;
initflag = MC(1).initflag;
if w(1,3) ~= 1
    vsflag = 1;
else
    vsflag = 0;
end
nlay = MC(1).nlay;
ray = MC(1).ray;
segmax = length(ray(end).path);
numph = length(ray);
phaselist = zeros(segmax,2,numph);
for i = 1:length(ray)
    nseg(i) = length(ray(i).path);
    phaselist(1:nseg(i),1,i) = ray(i).path;
    phaselist(1:nseg(i),2,i) = ray(i).type;
end
if nlay > 6
    nseg2 = nseg;
    phaselist2 = phaselist;
else
filename = ['raylist_' num2str(nlay) '_layers.txt'];
[ray] = readphaselist(filename);
segmax = length(ray(end).path);
numph = length(ray);
phaselist2 = zeros(segmax,2,numph);
for i = 1:length(ray)
    nseg2(i) = length(ray(i).path);
    phaselist2(1:nseg2(i),1,i) = ray(i).path;
    phaselist2(1:nseg2(i),2,i) = ray(i).type;
end
end

% Load data for inversion
% -----------------------
load(input_file)
clear geom
time = linspace(0,(param(2)-1).*param(3),param(2))+param(6);
zz = -25:0.25:150;
hm(:,zz>param(9)) = [];
zz(zz>param(9)) = [];
hm(:,zz<param(8)) = [];
zz(zz<param(8)) = [];
bootm0 = hm'./max(abs(hm(:)));
% if length(param)>9
    rrf = zeros(length(RF),length(time));
    trf = rrf;
    for i = 1:length(RF)
        rrf(i,:) = interp1(RF(i).time,RF(i).rad,time)./max(abs(RF(i).rad));
        trf(i,:) = interp1(RF(i).time,RF(i).tra,time)./max(abs(RF(i).rad));
        geom(i,1) = mod(RF(i).baz,360);
        geom(i,2) = RF(i).slow;
    end
    rrf(isnan(rrf)) = 0;
    trf(isnan(trf)) = 0;
    rrf(:,time>param(10)) = [];
    trf(:,time>param(10)) = [];
    time(time>param(10)) = [];
%     geom3 = geom;
%     geom3(end+1:end+length(geom2),:) = geom2;
% else
%     geom3 = geom;
% end
% geom3 = geom;
% z_out = -25:0.25:150.25;
ldepth = [20, 35, 300]; % depths of interfaces
velp = [5.8, 6.5, 8]; % P wave velocity
vels = [3.4, 3.75, 4.5]; % S wave velocity
for i = 1:length(geom(:,1))
    for j = 1:length(zz)
        t_ref{i}(j) = z2t(zz(j),ldepth,velp,vels,geom(i,2));
        % if j == 1
        %     disp(geom3(i,2))
        % end
    end
end
% for j = 1:length(z_out)
%     t_ref{length(geom3(:,2))+1}(j) = z2t(z_out(j),ldepth,velp,vels,0.05);
% end

% Character vector for the parameters
parnames{1} = 'Thickness';
parnames{2} = 'P-wave velocity';
parnames{3} = 'S-wave velocity';
parnames{4} = 'Density';
parnames{5} = 'Anisotropy';
parnames{6} = 'Symmetry axis';
parnames{7} = 'Plunge';
parnames{8} = 'Strike';
parnames{9} = 'Dip';
parnames{10} = 'Velocity contrast';
% parunit{1} = 'm';
% parunit{2} = 'm/s';
% parunit{3} = 'm/s';
% parunit{4} = 'g/m^3';
parunit{1} = 'km';
parunit{2} = 'km/s';
parunit{3} = 'km/s';
parunit{4} = 'kg/m^3';
parunit{5} = 'pct';
parunit{6} = 'deg';
parunit{7} = 'deg';
parunit{8} = 'deg';
parunit{9} = 'deg';
parunit{10} = '';

% Evaluation of the results
% -------------------------
for i = 1:Nmc
    if i == 1
        if isfield(MC(i).post,'full')
        ll = MC(i).post.full.L1';
        res = MC(i).post.full.stat';
        hk = MC(i).post.full.Hk';
        err = MC(i).post.full.rms';
        else
        ll = MC(i).post.ll; 
        res = MC(i).post.stat(1:end-1,:); 
        end
    else
        if isfield(MC(i).post,'full')
        ll = [MC(i).post.full.L1' ll];
        res = [MC(i).post.full.stat' res];
        hk = [MC(i).post.full.Hk' hk];
        err = [MC(i).post.full.rms' err];
        else
        ll = [MC(i).post.ll ll]; 
        res = [MC(i).post.stat(1:end-1,:) res];
        end
    end
end
ll(isnan(ll)) = 0;
if isfield(MC(1).post,'full')
ioi = hk~=0;
hk(ioi) = hk(ioi)-min(hk(ioi));
hk(ioi) = hk(ioi)./max(hk(ioi));
hk(~ioi) = 0;
hk(isnan(hk)) = 0;
ioi = err~=0;
err(ioi) = err(ioi)-min(err(ioi));
err(ioi) = err(ioi)./max(abs(err(ioi)));
err(~ioi) = 1;
err(isnan(err)) = 1;
bla = hk.*(1-err);
else
    bla = ll;
end
figure
plot(ll,bla,'.');
ll = bla;

ii = randi(length(ll),max(round(length(ll)/100),3*NN2),1);
for j = 1:length(ll(ii))
    disp([num2str(j) ' of ' num2str(length(ll(ii)))])
    temp = makemodel(reshape(res(:,ii(j)),size(w)),b);
    % if initflag
    % fm = max(max(abs(b.amax)),max(abs(b.amin)))./((max(max(abs(b.amax)),max(abs(b.amin)))-10).*(temp(:,1)./1000).^-0.5+10);
    % fm(fm<1) = (temp(fm<1,1)./1000).^-2;
    % temp(:,6:7) = temp(:,6:7)./fm;
    fm = max(max(abs(b.amax)),max(abs(b.amin)))./((max(max(abs(b.amax)),max(abs(b.amin)))-10).*(temp(:,1)./1000).^-0.25+10);
    fm(fm<1) = (temp(fm<1,1)./1000).^-2;
    temp(:,6:7) = temp(:,6:7)./fm;
    % end
    n = 1;
    for i = 1:length(temp(:,1))
        if i == 1
            hrt(j,n) = 0;
        else
            hrt(j,n) = hrt(j,n-1);
        end
        hrt(j,n+1) = hrt(j,n)+temp(i,1)./1000;
        htt(j,n) = temp(i,1)./1000;
        htt(j,n+1) = temp(i,1)./1000;
        rhort(j,n) = temp(i,2)./1000;
        rhort(j,n+1) = temp(i,2)./1000;
        vprt(j,n) = temp(i,3)./1000;
        vprt(j,n+1) = temp(i,3)./1000;
        vsrt(j,n) = temp(i,4)./1000;
        vsrt(j,n+1) = temp(i,4)./1000;
        art(j,n) = temp(i,6);
        art(j,n+1) = temp(i,6);
        phirt(j,n) = mod(temp(i,8),360);
        phirt(j,n+1) = mod(temp(i,8),360);
        plrt(j,n) = temp(i,9);
        plrt(j,n+1) = temp(i,9);
        srrt(j,n) = mod(temp(i,10),360);
        srrt(j,n+1) = mod(temp(i,10),360);
        dprt(j,n) = temp(i,11);
        dprt(j,n+1) = temp(i,11);
        n = n+2;
    end
end
hrt(:,end) = hrt(:,end).*2;
phirt(plrt<0) = phirt(plrt<0)+180;
plrt(plrt<0) = -plrt(plrt<0);


% keyboard
[~,ii] = sort(ll,'descend');
Nres = min(NN1,round(length(ii)./10.*9)); % Number of results accepted from all trials (Nmc*N)
I = ii(1:Nres);
loi = ll(I);
if isfield(MC(1).post,'full')
hkoi = hk(I);
hkoi = hkoi-min(hkoi);
hkoi = hkoi./max(abs(hkoi));
erroi = err(I);
erroi = erroi-min(erroi);
erroi = erroi./max(abs(erroi));
bloi = hkoi.*(1-erroi);
else
    bloi = loi;
end
% bloi = bla(I);
figure
plot(loi,bloi,'.');

for j = 1:length(ll(I))
    res2(:,:,j) = makemodel(reshape(res(:,I(j)),size(w)),b);
    % if initflag
    % fm = max(max(abs(b.amax)),max(abs(b.amin)))./((max(max(abs(b.amax)),max(abs(b.amin)))-10).*(res2(:,1,j)./1000).^-0.5+10);
    % fm(fm<1) = (res2(fm<1,1,j)./1000).^-2;
    % res2(:,6:7,j) = res2(:,6:7,j)./fm;
    temp = res2(:,:,j);
    fm = max(max(abs(b.amax)),max(abs(b.amin)))./((max(max(abs(b.amax)),max(abs(b.amin)))-10).*(temp(:,1)./1000).^-0.25+10);
    fm(fm<1) = (temp(fm<1,1)./1000).^-2;
    res2(:,6:7,j) = res2(:,6:7,j)./fm;
    % end
    % grvec(:,j) = res(:,I(j));
    temp2(:,1) = temp(:,1);
    temp2(:,2) = temp(:,4);
    temp2(:,3) = temp(:,6)>0;
    temp2(:,4) = temp(:,6)==0;
    temp2(:,5) = abs(temp(:,6));
    temp2(:,6) = temp(:,9);
    temp2(:,7) = temp(:,11);
    temp2(:,8) = ll(I(j));
    temp2(:,9) = sind(temp(:,8));
    temp2(:,10) = cosd(temp(:,8));
    temp2(:,11) = sind(temp(:,10));
    temp2(:,12) = cosd(temp(:,10));
    
    if j == 1
        gr2norm = zeros(size(temp2));
        gr2norm(:,1:8) = 1;
        gr2norm = gr2norm(:);
    end
    grvec(:,j) = temp2(:);
    n = 1;
    for i = 1:length(res2(:,1,j))
        disp([num2str(j) ' of ' num2str(length(ll(I))) ' - ' num2str(i) ' of ' num2str(length(res2(:,1,j)))]);
        if i == 1
            hr(j,n) = 0;
        else
            hr(j,n) = hr(j,n-1)+0.0001;
        end
        hr(j,n+1) = hr(j,n)+res2(i,1,j)./1000;
        ht(j,n) = res2(i,1,j)./1000;
        ht(j,n+1) = res2(i,1,j)./1000;
        rhor(j,n) = res2(i,2,j)./1000;
        rhor(j,n+1) = res2(i,2,j)./1000.*1.0001;
        vpr(j,n) = res2(i,3,j)./1000;
        vpr(j,n+1) = res2(i,3,j)./1000.*1.0001;
        vsr(j,n) = res2(i,4,j)./1000;
        vsr(j,n+1) = res2(i,4,j)./1000.*1.0001;
        ar(j,n) = res2(i,6,j);
        ar(j,n+1) = res2(i,6,j).*1.0001;
        phir(j,n) = mod(res2(i,8,j),360);
        phir(j,n+1) = mod(res2(i,8,j),360).*1.0001;
        plr(j,n) = res2(i,9,j);
        plr(j,n+1) = res2(i,9,j).*1.0001;
        srr(j,n) = mod(res2(i,10,j),360);
        srr(j,n+1) = mod(res2(i,10,j),360).*1.0001;
        dpr(j,n) = res2(i,11,j);
        dpr(j,n+1) = res2(i,11,j).*1.0001;
        n = n+2;
    end
end
hr(:,end) = hr(:,end).*2;
phir(plr<0) = phir(plr<0)+180;
plr(plr<0) = -plr(plr<0);

% res0 = [6400 2500 4200 2400 1 0 0 0 0 0 0;
%     4000 2700 5600 3200 0 -10 -10 45 20 90 10;
%     25000 2900 6300 3700 0 -5 -5 0 0 90 10;
%     0 3300 8000 4700 1 0 0 0 0 90 10];

% temp = geom3(:,2);
% temp(1./temp<9) = 1./9.;
% geom3(:,2) = temp;
% [hdc,rr,tt,wr,wt] = forwardHDray4(res2(:,:,loi==max(loi)),param,geom3,phaselist2,nseg2);
% [hdc,rr,tt,wr,wt] = forwardHDray4(res2(:,:,hkoi==max(hkoi)),param,geom3,phaselist2,nseg2);
% [hdc,rr,tt,wr,wt] = forwardHDray4(res0,param,geom3,phaselist2,nseg2);

% [hdc,rr,tt,wr,wt] = forwardHDray2(res2(:,:,bloi==max(bloi)),param,geom,t_ref,phaselist2,nseg2);
if ~isempty(modpath)
    fid = fopen(modpath,'rt');
    C = textscan(fid,'%f %f %f %f %f %f %f %f %f %f','Delimiter',',');
    fclose(fid);
    mho = C{1,1}./1000;
    mvpo = C{1,3}./1000;
    mvso = C{1,4}./1000;
    msro = C{1,5};
    mdpo = C{1,6};
    mao = C{1,8};
    mphio = C{1,9};
    mplo = C{1,10};
    msro(1) = [];
    mdpo(1) = [];
    msro(end+1) = 0;
    mdpo(end+1) = 0;
    mphio(mplo<0) = mphio(mplo<0)+180;
    mplo(mplo<0) = -mplo(mplo<0);
    n = 1;
    for i = 1:length(mho)
        if i == 1
            ho(n) = 0;
        else
            ho(n) = ho(n-1)+0.0001;
        end
        ho(n+1) = ho(n)+mho(i);
        vpo(n) = mvpo(i);
        vpo(n+1) = mvpo(i).*1.0001;
        vso(n) = mvso(i);
        vso(n+1) = mvso(i).*1.0001;
        ao(n) = mao(i);
        ao(n+1) = mao(i).*1.0001;
        n = n+2;
    end
    ho(:,end) = ho(:,end).*2;
end
if ~isempty(modpath)
    N = length(mho);
    % modelo = res2(:,:,1);
    modelo = zeros(N,11);
    modelo(:,1) = mho*1000;
    modelo(:,2) = C{1,2};
    modelo(:,3) = C{1,3};
    modelo(:,4) = C{1,4};
    modelo(:,5) = mao==0;
    modelo(:,6) = mao;
    modelo(:,7) = mao;
    modelo(:,8) = mphio;
    modelo(:,9) = mplo;
    modelo(2:end,10) = msro(1:end-1);
    modelo(2:end,11) = mdpo(1:end-1);
    
end
% keyboard

% figure; subplot(1,3,1); imagesc(bootm0(1:500,1:5)'); caxis([-1 1]); subplot(1,3,2); imagesc(hdc(1:500,1:5)'); caxis([-1 1]); subplot(1,3,3); imagesc(hdc(1:500,1:5)'-bootm0(1:500,1:5)'); caxis([-1 1])

% if exist('RF','var')
% M2 = size(wr);
% figure; 
% subplot(1,2,1);
% hold on
% for i = 1:length(RF)
%     plot(time,rrf(i,:)./max(abs(rrf(i,:)))+(i-1),'k','LineWidth',1.1)
%     plot(time,rr(i,:)./max(abs(rr(i,:)))+(i-1),'g')
%     plot(time,wr(i,:)+(i-1),'r--')
% end
% axis([-2 25 -1 length(RF)+1])
% subplot(1,2,2);
% hold on
% for i = 1:length(RF)
%     plot(time,trf(i,:)./max(abs(rrf(i,:)))+(i-1),'k','LineWidth',1.1)
%     plot(time,tt(i,:)./max(abs(rr(i,:)))+(i-1),'g')
%     plot(time,wt(i,:)+(i-1),'r--')
% end
% axis([-2 25 -1 length(RF)+1])
% % keyboard
% end

N1 = size(grvec,1);
grveco = grvec;
n = 1;
for i = 1:N1
    if std(grveco(i,:))<eps*1000
        todel(n) = i;
        n = n+1;
        continue
    else
        if gr2norm(i) == 1
            temp = grvec(i,:);
            temp = temp-min(temp(:));
            temp = temp./max(temp(:));
            grvec(i,:) = temp;
        end
    end
end
if n > 1
    grvec(todel,:) = [];
end
N2 = size(grvec,2);
coi = zeros(N2);
for i = 1:N2
%     for j = 1:N2
        coi(i,:) = sqrt(sum((grvec(:,i)-grvec(:,:)).^2,1));
%     end
    disp(num2str(i))
end
coi = coi-min(coi(:));
coi = coi./max(coi(:));
dismat = coi;
% Z = linkage(dismat,'complete');
Z = linkage(dismat,'ward');
lt = 0; ng = Nres; Ngroups = 3;
while ng(lt==max(lt)) > min(NN2,length(bloi)/10)%min(1000,length(ll)/15)
Ngroups = Ngroups+1;
groups = cluster(Z,'MaxClust',Ngroups);
indx = 1:length(groups);
groupsut = unique(groups);
jj = [];
for i = 1:length(groupsut)
    jj = [jj indx(groupsut(i)==groups)];
    lt(i) = max(bloi(groupsut(i)==groups));
    ng(i) = sum(groups==groupsut(i));
end
end
Ngroups = Ngroups-1;
groups = cluster(Z,'MaxClust',Ngroups);
indx = 1:length(groups);
groupsut = unique(groups);
jj = [];
for i = 1:length(groupsut)
    jj = [jj indx(groupsut(i)==groups)];
    lt(i) = max(bloi(groupsut(i)==groups));
    ng(i) = sum(groups==groupsut(i));
end

% keyboard
% for i = 1:length(groupsut)
%     ng(i) = sum(groups==groupsut(i));
% end
% [~,ii] = max(ng);
if ~isscalar(groupsut)
[~,ii] = max(lt);
i2del = indx(groups~=groupsut(ii));

I(i2del) = [];
hr(i2del,:) = [];
ht(i2del,:) = [];
rhor(i2del,:) = [];
vpr(i2del,:) = [];
vsr(i2del,:) = [];
ar(i2del,:) = [];
phir(i2del,:) = [];
plr(i2del,:) = [];
srr(i2del,:) = [];
dpr(i2del,:) = [];

bloi(i2del) = [];
loi(i2del) = [];
end

alph = bloi;
alph = alph-min(alph);
alph = alph./max(alph)/length(I)/2;
alph = alph.^0.5;
% alph = alph.^0.2;

% Make probability density of vp
% ------------------------------
hh = 0:0.05:(max(hr(:,end-1))+10);
pp = (min(vpr(:))-0.5):0.02:(max(vpr(:))+0.5);
ppi = 1:length(pp);
HP = zeros(length(hh),length(pp));
for j = 1:length(loi)
    vpi = round(interp1(pp,ppi,vpr(j,:)));
    vpi(isnan(vpi)) = ppi(end);
    vpi2 = round(interp1(hr(j,:),vpi,hh,'nearest'));
    vpi2(isnan(vpi2)) = ppi(end);
    try
        pos = sub2ind(size(HP),1:length(hh),vpi2);
    catch
        keyboard
    end
    HP(pos) = HP(pos)+1;
end
HP(HP>max(HP(round(length(hh)/2),:))) = max(HP(round(length(hh)/2),:));
Nc = 64;
amp = linspace(min(HP(HP>0)),max(HP(:)),Nc+1);
indx = 1:length(HP(:));
colp = hot(Nc+5);
for j = 1:Nc
    pos = indx(HP>amp(j)&HP<=amp(j+1));
    [cop(j).I,cop(j).J] = ind2sub(size(HP),pos);
    cop(j).h = hh(cop(j).I);
    cop(j).p = pp(cop(j).J);
end
% Make probability density of a
% -----------------------------
% vpar = vpr.*(1+ar./100);
vpar = ar/10;
% vpaa = (min(vpar(:))-0.5):0.02:(max(vpar(:))+0.5);
vpaa = linspace(min(vpar(:))-0.5,max(vpar(:))+0.5,500);
vpar(abs(ar)==0) = nan;
vpai = 1:length(vpaa);
HA = zeros(length(hh),length(vpaa));
for j = 1:length(loi)
    vpari = round(interp1(vpaa,vpai,vpar(j,:)));
    vpari2 = round(interp1(hr(j,:),vpari,hh,'nearest'));
    hi = 1:length(hh);
    pos = sub2ind(size(HA),hi(~isnan(vpari2)),vpari2(~isnan(vpari2)));
    HA(pos) = HA(pos)+1;
end
% HA(HA>max(HA(round(length(hh)/2),:))) = max(HA(round(length(hh)/2),:));
Nc = 64;
ioi = HA==0;
if sum(~ioi(:))==0
    amp = linspace(-20,20,Nc+1);
else
    amp = linspace(min(HA(~ioi)),max(HA(~ioi)),Nc+1);
end
indx = 1:length(HA(:));
% cola = cool(Nc);
cola = sky(Nc);
for j = 1:Nc
    pos = indx(HA>amp(j)&HA<=amp(j+1));
    [coa(j).I,coa(j).J] = ind2sub(size(HA),pos);
    coa(j).h = hh(coa(j).I);
    coa(j).a = vpaa(coa(j).J);
end
% Make probability density of vs
% ------------------------------
ss = (min(vsr(:))-0.5):0.02:(max(vsr(:))+0.5);
si = 1:length(ss);
HS = zeros(length(hh),length(ss));
for j = 1:length(loi)
    vsi = round(interp1(ss,si,vsr(j,:)));
    vsi(isnan(vsi)) = si(end);
    vsi2 = round(interp1(hr(j,:),vsi,hh,'nearest'));
    vsi2(isnan(vsi2)) = vsi(end);
    pos = sub2ind(size(HS),1:length(hh),vsi2);
    HS(pos) = HS(pos)+1;
end
HS(HS>max(HS(round(length(hh)/2),:))) = max(HS(round(length(hh)/2),:));
Nc = 64;
amp = linspace(min(HS(HS>0)),max(HS(:)),Nc+1);
indx = 1:length(HS(:));
% cols = bone(Nc+5);
cols = parula(Nc);
for j = 1:Nc
    pos = indx(HS>amp(j)&HS<=amp(j+1));
    [coss(j).I,coss(j).J] = ind2sub(size(HS),pos);
    coss(j).h = hh(coss(j).I);
    coss(j).s = ss(coss(j).J);
end

% Prepare statistics for parameters
for i = 1:nlay-1
    stat(i).vp = vpr(:,2*i-1);
    stat(i).mvp = vpr(1,2*i-1);%mean(vpr(:,2*i-1));
    stat(i).mvp2 = mean(vpr(:,2*i-1));
    stat(i).dvp = 2*sqrt(var(vpr(:,2*i-1)));
    stat(i).nvp = length(vpr(:,2*i-1));
    stat(i).h = ht(:,2*i-1);
    stat(i).mh = ht(1,2*i-1);%mean(ht(:,2*i-1));
    stat(i).mh2 = mean(ht(:,2*i-1));
    stat(i).dh = 2*sqrt(var(ht(:,2*i-1)));
    stat(i).nh = length(ht(:,2*i-1));
    stat(i).vs = vsr(:,2*i-1);
    stat(i).mvs = vsr(1,2*i-1);%mean(vsr(:,2*i-1));
    stat(i).mvs2 = mean(vsr(:,2*i-1));
    stat(i).dvs = 2*sqrt(var(vsr(:,2*i-1)));
    stat(i).nvs = length(vsr(:,2*i-1));
    stat(i).rho = rhor(:,2*i-1);
    stat(i).mrho = rhor(1,2*i-1);%mean(rhor(:,2*i-1));
    stat(i).mrho2 = mean(rhor(:,2*i-1));
    stat(i).drho = 2*sqrt(var(rhor(:,2*i-1)));
    stat(i).nrho = length(rhor(:,2*i-1));
    stat(i).a = ar(:,2*i-1);
    stat(i).ma = ar(1,2*i-1);%mean(ar(:,2*i-1));
    stat(i).ma2 = mean(ar(:,2*i-1));
    stat(i).da = 2*sqrt(var(ar(:,2*i-1)));
    stat(i).na = length(ar(:,2*i-1));
    stat(i).phi = mod(phir(:,2*i-1),360);
    stat(i).nphi = length(phir(:,2*i-1));
    vartemp = mod(phir(:,2*i-1),360);
    vartemp2 = vartemp;
    vartemp2(vartemp2>180) = vartemp2(vartemp2>180)-360;
    dvartemp1 = 2*sqrt(var(vartemp));
    dvartemp2 = 2*sqrt(var(vartemp2));
    if dvartemp1<dvartemp2
        stat(i).mphi = vartemp(1);%mean(vartemp);
        stat(i).mphi2 = mean(vartemp);
        stat(i).dphi = dvartemp1;
    else
        stat(i).mphi = mod(vartemp2(1),360);%mod(mean(vartemp2),180);
        stat(i).mphi2 = mod(mean(vartemp2),360);
        stat(i).dphi = dvartemp2;
    end
    stat(i).pl = plr(:,2*i-1);
    stat(i).mpl = plr(1,2*i-1);%mean(plr(:,2*i-1));
    stat(i).mpl2 = mean(plr(:,2*i-1));
    stat(i).dpl = 2*sqrt(var(plr(:,2*i-1)));
    stat(i).npl = length(plr(:,2*i-1));
    stat(i).nsr = length(srr(:,2*(i+1)-1));
    stat(i).sr = mod(srr(:,2*(i+1)-1),360);
    vartemp = mod(srr(:,2*(i+1)-1),360);
    vartemp2 = vartemp;
    vartemp2(vartemp2>180) = vartemp2(vartemp2>180)-360;
    dvartemp1 = 2*sqrt(var(vartemp));
    dvartemp2 = 2*sqrt(var(vartemp2));
    if dvartemp1<dvartemp2
        stat(i).msr = vartemp(1);%mean(vartemp);
        stat(i).msr2 = mean(vartemp);
        stat(i).dsr = dvartemp1;
    else
        stat(i).msr = mod(vartemp2(1),360);%mod(mean(vartemp2),360);
        stat(i).msr2 = mod(mean(vartemp2),360);
        stat(i).dsr = dvartemp2;
    end
    stat(i).dp = dpr(:,2*(i+1)-1);
    stat(i).mdp = dpr(1,2*(i+1)-1);%mean(dpr(:,2*(i+1)-1));
    stat(i).mdp2 = mean(dpr(:,2*(i+1)-1));
    stat(i).ddp = 2*sqrt(var(dpr(:,2*(i+1)-1)));
    stat(i).ndp = length(dpr(:,2*(i+1)-1));
end
stat(nlay).vp = vpr(:,2*nlay-1);
stat(nlay).mvp = vpr(1,2*nlay-1);%mean(vpr(:,2*nlay-1));
stat(nlay).mvp2 = mean(vpr(:,2*nlay-1));
stat(nlay).dvp = 2*sqrt(var(vpr(:,2*nlay-1)));
stat(nlay).nvp = length(vpr(:,2*nlay-1));
stat(nlay).h = ht(:,2*nlay-1);
stat(nlay).mh = ht(1,2*nlay-1);%mean(ht(:,2*nlay-1));
stat(nlay).mh2 = mean(ht(:,2*nlay-1));
stat(nlay).dh = 2*sqrt(var(ht(:,2*nlay-1)));
stat(nlay).nh = length(ht(:,2*nlay-1));
stat(nlay).vs = vsr(:,2*nlay-1);
stat(nlay).mvs = vsr(1,2*nlay-1);%mean(vsr(:,2*nlay-1));
stat(nlay).mvs2 = mean(vsr(:,2*nlay-1));
stat(nlay).dvs = 2*sqrt(var(vsr(:,2*nlay-1)));
stat(nlay).nvs = length(vsr(:,2*nlay-1));
stat(nlay).rho = rhor(:,2*nlay-1);
stat(nlay).mrho = rhor(1,2*nlay-1);%mean(rhor(:,2*nlay-1));
stat(nlay).mrho2 = mean(rhor(:,2*nlay-1));
stat(nlay).drho = 2*sqrt(var(rhor(:,2*nlay-1)));
stat(nlay).nrho = length(rhor(:,2*nlay-1));
stat(nlay).a = ar(:,2*nlay-1);
stat(nlay).ma = 0;
stat(nlay).ma2 = 0;
stat(nlay).da = 0;
stat(nlay).na = length(ar(:,2*nlay-1));
stat(nlay).phi = mod(phir(:,2*nlay-1),360);
stat(nlay).nphi = length(phir(:,2*nlay-1));
stat(nlay).mphi = 0;
stat(nlay).mphi2 = 0;
stat(nlay).dphi = 0;
stat(nlay).pl = plr(:,2*nlay-1);
stat(nlay).mpl = 0;
stat(nlay).mpl2 = 0;
stat(nlay).dpl = 0;
stat(nlay).npl = length(plr(:,2*nlay-1));
stat(nlay).nsr = 0;
stat(nlay).sr = 0;
stat(nlay).msr = 0;
stat(nlay).msr2 = 0;
stat(nlay).dsr = 0;
stat(nlay).dp = 0;
stat(nlay).mdp = 0;
stat(nlay).mdp2 = 0;
stat(nlay).ddp = 0;
stat(nlay).ndp = 0;

if ~isempty(modpath)
    fid = fopen(modpath,'rt');
    C = textscan(fid,'%f %f %f %f %f %f %f %f %f %f','Delimiter',',');
    fclose(fid);
    mho = C{1,1}./1000;
    mvpo = C{1,3}./1000;
    mvso = C{1,4}./1000;
    msro = C{1,5};
    mdpo = C{1,6};
    mao = C{1,8};
    mphio = C{1,9};
    mplo = C{1,10};
    msro(1) = [];
    mdpo(1) = [];
    msro(end+1) = 0;
    mdpo(end+1) = 0;
    mphio(mplo<0) = mphio(mplo<0)+180;
    mplo(mplo<0) = -mplo(mplo<0);
    n = 1;
    for i = 1:length(mho)
        if i == 1
            ho(n) = 0;
        else
            ho(n) = ho(n-1)+0.0001;
        end
        ho(n+1) = ho(n)+mho(i);
        vpo(n) = mvpo(i);
        vpo(n+1) = mvpo(i).*1.0001;
        vso(n) = mvso(i);
        vso(n+1) = mvso(i).*1.0001;
        ao(n) = mao(i);
        ao(n+1) = mao(i).*1.0001;
        n = n+2;
    end
    ho(:,end) = ho(:,end).*2;
end

n = 1;
for i = 1:length(stat)
    if i == 1
        hf(n) = 0;
        hf2(n) = 0;
    else
        hf(n) = hf(n-1)+0.0001;
        hf2(n) = hf2(n-1)+0.0001;
    end
    hf(n+1) = hf(n)+stat(i).mh;
    hf2(n+1) = hf2(n)+stat(i).mh2;
%     rhof(n) = res2(i,2,j)./1000;
%     rhof(n+1) = res2(i,2,j)./1000.*1.0001;
    vpf(n) = stat(i).mvp;
    vpf(n+1) = stat(i).mvp.*1.0001;
    vsf(n) = stat(i).mvs;
    vsf(n+1) = stat(i).mvs.*1.0001;
    af(n) = stat(i).ma;
    af(n+1) = stat(i).ma.*1.0001;
    vpf2(n) = stat(i).mvp2;
    vpf2(n+1) = stat(i).mvp2.*1.0001;
    vsf2(n) = stat(i).mvs2;
    vsf2(n+1) = stat(i).mvs2.*1.0001;
    af2(n) = stat(i).ma2;
    af2(n+1) = stat(i).ma2.*1.0001;
%     phif(n) = mod(res2(i,8,j),180);
%     phif(n+1) = mod(res2(i,8,j),180).*1.0001;
%     plf(n) = res2(i,9,j);
%     plf(n+1) = res2(i,9,j).*1.0001;
%     srf(n) = mod(res2(i,10,j),360);
%     srf(n+1) = mod(res2(i,10,j),360).*1.0001;
%     dpf(n) = res2(i,11,j);
%     dpf(n+1) = res2(i,11,j).*1.0001;
    n = n+2;
end
hf(:,end) = hf(:,end).*2;
hf2(:,end) = hf2(:,end).*2;


FS = 9;
fig = figure('units','normalized','outerposition',[0 0 1 1]);
% subplot(nlay+2,6,sort([2:4:(6*nlay+2) 3:4:(6*nlay+2)]))
ax0 = subplot(nlay+1,8,[3:8:((nlay+1)*8)' 4:8:((nlay+1)*8)' 5:8:((nlay+1)*8)' 6:8:((nlay+1)*8)']);
hold on
if ~vsflag
plot(vprt',hrt','Color',[0.7 0.7 0.7])
end
plot(art'./10,hrt','Color',[0.7 0.7 0.7])
plot(vsrt',hrt','Color',[0.7 0.7 0.7])
set(gca, 'YDir','reverse')
for j = 1:Nc
    if ~vsflag
    scatter(cop(j).p,cop(j).h,2,colp(j,:),'filled')
    end
    scatter(coa(j).a,coa(j).h,2,cola(j,:),'filled')
    scatter(coss(j).s,coss(j).h,2,cols(j,:),'filled')
end
for i = 1:nlay
plot([min([min(vprt(:))-0.5 min(vpar(:))-1]) max(vprt(:))+0.5],[hf(i*2) hf(i*2)],'Color',[80 220 100]./255,'LineWidth',0.5)
end
if ~isempty(modpath)
    if ~vsflag
        % plot(vpo,ho,'-','Color',([0 100 255])./255,'LineWidth',1.5)
        plot(vpo,ho,'-','Color',([125 100 0])./255,'LineWidth',1.5)
    end
% plot(ao./10,ho,'-','Color',([0 100 255])./255,'LineWidth',1.5)
% plot(vso,ho,'-','Color',([0 100 255])./255,'LineWidth',1.5)
plot(ao./10,ho,'-','Color',([125 100 0])./255,'LineWidth',1.5)
plot(vso,ho,'-','Color',([125 100 0])./255,'LineWidth',1.5)
end
if ~bfflag
    if ~vsflag
        % plot(vpf2,hf2,'--','Color',([80 255 100])./255,'LineWidth',1.5)
        plot(vpf2,hf2,'--','Color',([255 80 100])./255,'LineWidth',1.5)
    end
% plot(af2./10,hf2,'--','Color',([80 255 100])./255,'LineWidth',1.5)
% plot(vsf2,hf2,'--','Color',([80 255 100])./255,'LineWidth',1.5)
plot(af2./10,hf2,'--','Color',([255 80 100])./255,'LineWidth',1.5)
plot(vsf2,hf2,'--','Color',([255 80 100])./255,'LineWidth',1.5)
else
    if ~vsflag
        % plot(vpf,hf,'--','Color',[80 255 100]./255,'LineWidth',1.5)
        plot(vpf,hf,'--','Color',[255 80 100]./255,'LineWidth',1.5)
    end
% plot(af./10,hf,'--','Color',[80 255 100]./255,'LineWidth',1.5)
% plot(vsf,hf,'--','Color',[80 255 100]./255,'LineWidth',1.5)
plot(af./10,hf,'--','Color',[255 80 100]./255,'LineWidth',1.5)
plot(vsf,hf,'--','Color',[255 80 100]./255,'LineWidth',1.5)
end
text(0,-1,'Anisotropy in [%]','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',FS)
if vsflag
    text(vsf(1)-0.01,-1,'S-velocity in [km/s]','HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',FS)
    plot([min([min(vsrt(:))-0.5 min(vpar(:))-1]) max(vsrt(:))+0.5],[0 0],'k','LineWidth',1)
else
text(vsf(1)-0.01,-1,'P- and S-velocity in [km/s]','HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',FS)
plot([min([min(vprt(:))-0.5 min(vpar(:))-1]) max(vprt(:))+0.5],[0 0],'k','LineWidth',1)
end
plot([0 0],[0 max(hr(:,end-1))+5],'k','LineWidth',1)
if vsflag
% axis([min([min(vsrt(:))-0.5 min(vpar(:))-1]) max(vsrt(:))+0.5 -5 max(hr(:,end-1))+5])
axis([-3 max(vsrt(:))+0.5 -5 max(hr(:,end-1))+5])
else
% axis([min([min(vprt(:))-0.5 min(vpar(:))-1]) max(vprt(:))+0.5 -5 max(hr(:,end-1))+5])
axis([-3 max(vprt(:))+0.5 -5 max(hr(:,end-1))+5])
end
% ylh = ylabel('Depth in [km]','FontSize',FS);
% ylh.Position(2) = ylh.Position(2) + 0.3;
xl = xticklabels(gca);
for i = 1:length(xl)
    temp = str2double(xl{i});
    if temp < 2.5
        xl{i} = ['\color{blue} ' num2str(temp*10) ' \color{black}/' num2str(temp)];
    end
end
xticklabels(gca,xl);
yl = yticklabels(gca);
for i = 1:length(yl)
    temp = str2double(yl{i});
    if temp < 0
        yl{i} = '';
    end
end
yticklabels(gca,yl);
box on
grid on
set(gca,'FontSize',FS)
hold off

if exist('coa','var')
pos = get(ax0,'Position');
axcb1 = axes('Position',[pos(1)+0.01 pos(2)+pos(4) pos(3)/3-0.04 0.01]);
cb = colorbar(axcb1,'northoutside','AxisLocation','out');
colormap(axcb1,sky(Nc));
ylabel(cb,'Probability for a')
caxis([0 1])
set(axcb1,'visible','off')
end
if exist('coss','var')
pos = get(ax0,'Position');
axcb2 = axes('Position',[pos(1)+0.01+pos(3)/3 pos(2)+pos(4) pos(3)/3-0.04 0.01]);
cb = colorbar(axcb2,'northoutside','AxisLocation','out');
colormap(axcb2,parula(Nc));
ylabel(cb,'Probability for vs')
caxis([0 1])
set(axcb2,'visible','off')
end
if exist('cop','var')
pos = get(ax0,'Position');
axcb3 = axes('Position',[pos(1)+0.01+pos(3)/3*2 pos(2)+pos(4) pos(3)/3-0.04 0.01]);
cb = colorbar(axcb3,'northoutside','AxisLocation','out');
colormap(axcb3,hot(Nc));
ylabel(cb,'Probability for vp')
caxis([0 1])
set(axcb3,'visible','off')
end

for i = 1:nlay-1
subplot(nlay-1,8,(i-1)*8+1)
theta0 = [deg2rad(phirt(:,i*2))];% deg2rad((phirt(:,i*2)+180))];
ph0 = polarhistogram(theta0,25,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',.3);
hold on
theta = [deg2rad(phir(:,i*2))];% deg2rad((phir(:,i*2)+180))];
ph = polarhistogram(theta,25,'FaceColor','r','FaceAlpha',.3);
if ~isempty(modpath)
    polarplot(deg2rad([mphio(i) mphio(i)]),[0 max(ph.Values).*1.1],'Color',([0 100 255])./255,'LineWidth',1.5)
    % polarplot(deg2rad([mphio(i) mphio(i)]+180),[0 max(ph.Values).*1.1],'Color',([0 100 255])./255,'LineWidth',1.5)
end
if ~bfflag
polarplot(deg2rad([stat(i).mphi2 stat(i).mphi2]),[0 max(ph.Values).*1.1],'r','LineWidth',1.1)
% polarplot(deg2rad([stat(i).mphi2 stat(i).mphi2]+180),[0 max(ph.Values).*1.1],'r','LineWidth',1.1)
polarplot(deg2rad(stat(i).mphi2),max(ph.Values).*1.1,'r*','LineWidth',1.2,'MarkerSize',8)
else
polarplot(deg2rad([stat(i).mphi stat(i).mphi]),[0 max(ph.Values).*1.1],'r','LineWidth',1.1)
% polarplot(deg2rad([stat(i).mphi stat(i).mphi]+180),[0 max(ph.Values).*1.1],'r','LineWidth',1.1)
polarplot(deg2rad(stat(i).mphi),max(ph.Values).*1.1,'r*','LineWidth',1.2,'MarkerSize',8)
end
set(gca,'ThetaDir','clockwise','ThetaZeroLocation','Top','ThetaAxisUnits','degrees');
dtheta = 45;
thetaticks(0:dtheta:360-dtheta)
rlim([0 max(ph.Values).*1.1])
rticks(linspace(0,max(ph.Values).*1.1,5))
rticklabels('')
set(gca,'FontSize',FS)
% text(0,max(ph.Values)*1.6,'N','FontSize',10,'FontWeight','bold')
hold off
end

for i = 1:nlay-1
subplot(nlay-1,8,(i-1)*8+2)
temp = plrt(:,i*2);
temp(temp<0) = temp(temp<0)+180;
theta0 = deg2rad(temp);
ph0 = polarhistogram(theta0,25,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',.3);
hold on
temp = plr(:,i*2);
temp(temp<0) = temp(temp<0)+180;
theta = deg2rad(temp);
ph = polarhistogram(theta,15,'FaceColor','r','FaceAlpha',.3);
rlim([0 max(ph.Values).*1.1])
set(gca,'ThetaDir','clockwise','ThetaZeroLocation','right','ThetaAxisUnits','degrees');
dtheta = 30;
thetaticks(0:dtheta:180)
tticks = 0:dtheta:180;
for ii = 1:length(tticks)
    if tticks(ii) < 90
        tl{ii} = num2str(tticks(ii));
    elseif tticks(ii) == 90
        tl{ii} = '';
    else
        tl{ii} = num2str(tticks(ii)-180);
    end
end
set(gca,'ThetaTickLabel',tl)
rticks(linspace(0,max(ph.Values).*1.1,5))
rticklabels('')
if ~isempty(modpath)
    if mplo(i)>0
    polarplot(deg2rad([mplo(i) mplo(i)]),[0 max(ph.Values).*1.1],'Color',([0 100 255])./255,'LineWidth',1.5)
    else
    polarplot(deg2rad([mplo(i) mplo(i)]+180),[0 max(ph.Values).*1.1],'Color',([0 100 255])./255,'LineWidth',1.5)
    end
end
polarplot(deg2rad([0 0]),[0 max(ph.Values).*1.1],'r','LineWidth',1.1)
polarplot(deg2rad([180 180]),[0 max(ph.Values).*1.1],'r','LineWidth',1.1)
polarplot(0,max(ph.Values).*1.1,'r*','LineWidth',1.2,'MarkerSize',8)
if stat(i).mpl>0
    if ~bfflag
    polarplot(deg2rad([stat(i).mpl2 stat(i).mpl2]),[0 max(ph.Values).*1.1],'g','LineWidth',1.2)
    else
    polarplot(deg2rad([stat(i).mpl stat(i).mpl]),[0 max(ph.Values).*1.1],'g','LineWidth',1.2)
    end
else
    if ~bfflag
    polarplot(deg2rad([stat(i).mpl2 stat(i).mpl2]+180),[0 max(ph.Values).*1.1],'g','LineWidth',1.2)
    else
    polarplot(deg2rad([stat(i).mpl stat(i).mpl]+180),[0 max(ph.Values).*1.1],'g','LineWidth',1.2)
    end
end
set(gca,'ThetaLim',[0 180])
set(gca,'FontSize',FS)
hold off
end

for i = 1:nlay-1
subplot(nlay-1,8,(i-1)*8+7)
theta0 = deg2rad(srrt(:,(i+1)*2));
ph0 = polarhistogram(theta0,25,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',.3);
hold on
theta = deg2rad(srr(:,(i+1)*2));
ph = polarhistogram(theta,25,'FaceColor','r','FaceAlpha',.3);
if ~isempty(modpath)
    polarplot(deg2rad([msro(i) msro(i)]),[0 max(ph.Values).*1.1],'Color',([0 100 255])./255,'LineWidth',1.5)
end
if ~bfflag
polarplot(deg2rad([stat(i).msr2 stat(i).msr2]),[0 max(ph.Values).*1.1],'r','LineWidth',1.1)
polarplot(deg2rad([stat(i).msr2 stat(i).msr2]+90),[0 max(ph.Values).*1.1],'r--','LineWidth',1.1)
polarplot(deg2rad(stat(i).msr2+90),max(ph.Values).*1.1,'r*','LineWidth',1.2,'MarkerSize',8)
else
polarplot(deg2rad([stat(i).msr stat(i).msr]),[0 max(ph.Values).*1.1],'r','LineWidth',1.1)
polarplot(deg2rad([stat(i).msr stat(i).msr]+90),[0 max(ph.Values).*1.1],'r--','LineWidth',1.1)
polarplot(deg2rad(stat(i).msr+90),max(ph.Values).*1.1,'r*','LineWidth',1.2,'MarkerSize',8)
end
set(gca,'ThetaDir','clockwise','ThetaZeroLocation','Top','ThetaAxisUnits','degrees');
dtheta = 45;
thetaticks(0:dtheta:360-dtheta)
rlim([0 max(ph.Values).*1.1])
rticks(linspace(0,max(ph.Values).*1.1,5))
rticklabels('')
% text(0,max(ph.Values)*1.6,'N','FontSize',10,'FontWeight','bold')
set(gca,'FontSize',FS)
hold off
end

% for i = 1:nlay-1
% subplot(nlay-1,8,(i-1)*8+8)
% temp = dprt(:,(i+1)*2);
% % temp(temp<0) = temp(temp<0)+180;
% theta0 = deg2rad(temp);
% ph0 = polarhistogram(theta0,10,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',.3);
% hold on
% temp = dpr(:,(i+1)*2);
% % temp(temp<0) = temp(temp<0)+180;
% theta = deg2rad(temp);
% ph = polarhistogram(theta,5,'FaceColor','r','FaceAlpha',.3);
% rlim([0 max(ph.Values).*1.1])
% set(gca,'ThetaDir','clockwise','ThetaZeroLocation','right','ThetaAxisUnits','degrees');
% dtheta = 15;
% thetaticks(0:dtheta:45)
% rticks(linspace(0,max(ph.Values).*1.1,5))
% rticklabels('')
% if ~isempty(modpath)
%     polarplot(deg2rad([mdpo(i) mdpo(i)]),[0 max(ph.Values).*1.1],'Color',([0 100 255])./255,'LineWidth',1.5)
% end
% polarplot(deg2rad([0 0]),[0 max(ph.Values).*1.1],'r--','LineWidth',1.1)
% polarplot(0,max(ph.Values).*1.1,'r*','LineWidth',1.2,'MarkerSize',8)
% if ~bfflag
% polarplot(deg2rad([stat(i).mdp2 stat(i).mdp2]),[0 max(ph.Values).*1.1],'g','LineWidth',1.2)
% else
% polarplot(deg2rad([stat(i).mdp stat(i).mdp]),[0 max(ph.Values).*1.1],'g','LineWidth',1.2)
% end
% set(gca,'ThetaLim',[0 90])
% set(gca,'FontSize',FS)
% hold off
% end

for i = 1:nlay-1
subplot(nlay-1,8,(i-1)*8+8)
temp = dprt(:,(i+1)*2);
% temp(temp<0) = temp(temp<0)+180;
theta0 = deg2rad(temp);
ph0 = polarhistogram(theta0*3,15,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',.3);
hold on
temp = dpr(:,(i+1)*2);
% temp(temp<0) = temp(temp<0)+180;
theta = deg2rad(temp);
ph = polarhistogram(theta*3,10,'FaceColor','r','FaceAlpha',.3);
rlim([0 max(ph.Values).*1.1])
set(gca,'ThetaDir','clockwise','ThetaZeroLocation','right','ThetaAxisUnits','degrees');
dtheta = 30;
thetaticks(0:dtheta:90)
rticks(linspace(0,max(ph.Values).*1.1,5))
rticklabels('')
if ~isempty(modpath)
    polarplot(deg2rad([mdpo(i) mdpo(i)].*3),[0 max(ph.Values).*1.1],'Color',([0 100 255])./255,'LineWidth',1.5)
end
polarplot(deg2rad([0 0]),[0 max(ph.Values).*1.1],'r--','LineWidth',1.1)
polarplot(0,max(ph.Values).*1.1,'r*','LineWidth',1.2,'MarkerSize',8)
if ~bfflag
polarplot(deg2rad([stat(i).mdp2 stat(i).mdp2].*3),[0 max(ph.Values).*1.1],'g','LineWidth',1.2)
else
polarplot(deg2rad([stat(i).mdp stat(i).mdp].*3),[0 max(ph.Values).*1.1],'g','LineWidth',1.2)
end
set(gca,'ThetaLim',[0 90])
set(gca,'FontSize',FS)
thtks = 0:dtheta:90;
for it = 1:length(thtks)
    thtkslbs{it} = num2str(round(thtks(it)/3));
end
thetaticklabels(thtkslbs)
hold off
end

print(fig,'-r300','-djpeg','-painters',[savepath '/results_full' addname '.jpg']);
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,'-r0','-dpdf','-vector','-bestfit',[savepath '/results_full' addname '.pdf']);
%print(fig,'-r0','-depsc2','-vector',[savepath '/results_full' addname '.eps']);
close(fig);


for i = 1:length(stat)
    if ~bfflag
        modelr(i,1) = stat(i).mh2.*1000;
        modelr(i,2) = stat(i).mrho2.*1000;
        modelr(i,3) = stat(i).mvp2.*1000;
        modelr(i,4) = stat(i).mvs2.*1000;
        modelr(i,5) = stat(i).ma2 == 0;
        modelr(i,6) = stat(i).ma2;
        modelr(i,7) = stat(i).ma2;
        modelr(i,8) = stat(i).mphi2;
        modelr(i,9) = stat(i).mpl2;
        if i == 1
            modelr(i,10) = 0;
            modelr(i,11) = 0;
        else
            modelr(i,10) = stat(i-1).msr2;
            modelr(i,11) = stat(i-1).mdp2;
        end
    else
        modelr(i,1) = stat(i).mh.*1000;
        modelr(i,2) = stat(i).mrho.*1000;
        modelr(i,3) = stat(i).mvp.*1000;
        modelr(i,4) = stat(i).mvs.*1000;
        modelr(i,5) = stat(i).ma == 0;
        modelr(i,6) = stat(i).ma;
        modelr(i,7) = stat(i).ma;
        modelr(i,8) = stat(i).mphi;
        modelr(i,9) = stat(i).mpl;
        if i == 1
            modelr(i,10) = 0;
            modelr(i,11) = 0;
        else
            modelr(i,10) = stat(i-1).msr;
            modelr(i,11) = stat(i-1).mdp;
        end
    end
end
fid = fopen([savepath '/modelres_' addname '.txt'],'wt');
for i = 1:length(stat)
fprintf(fid,'%f %f %f %f %f %f %f %f %f %f\n',modelr(i,1),modelr(i,2),modelr(i,3),modelr(i,4),modelr(i,10),modelr(i,11),modelr(i,5),modelr(i,6),modelr(i,8),modelr(i,9));
end
fclose(fid);

if ~isempty(modpath)
    % modelo = modelr;
    N = length(mho);
    modelo = zeros(N,11);
    modelo(:,1) = mho*1000;
    modelo(:,2) = C{1,2};
    modelo(:,3) = C{1,3};
    modelo(:,4) = C{1,4};
    modelo(:,5) = mao==0;
    modelo(:,6) = mao;
    modelo(:,7) = mao;
    modelo(:,8) = mphio;
    modelo(:,9) = mplo;
    modelo(2:end,10) = msro(1:end-1);
    modelo(2:end,11) = mdpo(1:end-1);
    filename = ['raylist_' num2str(N) '_layers.txt'];
    [ray] = readphaselist(filename);
    segmax = length(ray(end).path);
    numph = length(ray);
    phaselist3 = zeros(segmax,2,numph);
    for i = 1:length(ray)
        nseg3(i) = length(ray(i).path);
        phaselist3(1:nseg3(i),1,i) = ray(i).path;
        phaselist3(1:nseg3(i),2,i) = ray(i).type;
    end
    [bootmo,rto,tto] = forwardHDray2(modelo,param,geom,t_ref,phaselist3,nseg3);
end
[bootmr,rtr,ttr] = forwardHDray2(modelr,param,geom,t_ref,phaselist2,nseg2);
% bootmr(time>param(8),:) = [];
if ~isempty(geompath)
    fid = fopen([geompath '/bazlist.csv'],'rt');
    C = textscan(fid,'%f');
    fclose(fid);
    bazn = C{1,1};
    fid = fopen([geompath '/slowlist.csv'],'rt');
    C = textscan(fid,'%f');
    fclose(fid);
    slown = C{1,1};
    geomn(1,:) = bazn;
    geomn(2,:) = slown;
    for i = 1:length(geomn(1,:))
        for j = 1:length(zz)
            t_refn{i}(j) = z2t(zz(j),ldepth,velp,vels,geomn(2,i));
            % if j == 1
            %     disp(geom3(i,2))
            % end
        end
    end
    [bootmr] = forwardHDray2(modelr,param,geomn',t_refn,phaselist2,nseg2);
    if ~isempty(modpath)
        [bootmo] = forwardHDray2(modelo,param,geomn',t_refn,phaselist2,nseg2);
    end
end

if exist('RF','var')
[hdc,rr,tt,wr,wt] = forwardHDray2(res2(:,:,bloi==max(bloi)),param,geom,t_ref,phaselist2,nseg2);
M2 = size(wr);
fig = figure; 
subplot(1,2,1);
hold on
for i = 1:length(RF)
    plot(time,rrf(i,:)./max(abs(rrf(i,:)))+(i-1),'k','LineWidth',1.2)
    plot(time,rtr(i+M2(1)-length(RF),:)./max(abs(rtr(i+M2(1)-length(RF),:)))+(i-1),'Color',[0 174 0 255/2]./255,'LineWidth',2)
    if ~isempty(modpath)
    plot(time,rto(i+M2(1)-length(RF),:)./max(abs(rto(i+M2(1)-length(RF),:)))+(i-1),'r--','LineWidth',1.1)
    end
end
axis([-2 param(10) -1 length(RF)])
xlabel('t in seconds')
ylabel('RF-bins')
title('Radial')
subplot(1,2,2);
hold on
for i = 1:length(RF)
    plot(time,trf(i,:)./max(abs(rrf(i,:)))+(i-1),'k','LineWidth',1.2)
    plot(time,ttr(i+M2(1)-length(RF),:)./max(abs(rtr(i+M2(1)-length(RF),:)))+(i-1),'Color',[0 174 0 255/2]./255,'LineWidth',2)
    if ~isempty(modpath)
    plot(time,tto(i+M2(1)-length(RF),:)./max(abs(rto(i+M2(1)-length(RF),:)))+(i-1),'r--','LineWidth',1.1)
    end
end
axis([-2 param(10) -1 length(RF)])
xlabel('t in seconds')
title('Transverse')
print(fig,'-r300','-djpeg','-painters',[savepath '/RFbinsFinal' addname '.jpg']);
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,'-r0','-dpdf','-vector','-bestfit',[savepath '/RFbinsFinal' addname '.pdf']);
%print(fig,'-r0','-depsc2','-vector',[savepath '/RFbinsFinal' addname '.eps']);
% keyboard
end




% for j = 1:min(NN2,length(stat(1).vp))
%     clear modelt
% for i = 1:length(stat)
%     modelt(i,1) = stat(i).h(j).*1000;
%     modelt(i,2) = stat(i).rho(j).*1000;
%     modelt(i,3) = stat(i).vp(j).*1000;
%     modelt(i,4) = stat(i).vs(j).*1000;
%     modelt(i,5) = stat(i).a(j) == 0;
%     modelt(i,6) = stat(i).a(j);
%     modelt(i,7) = stat(i).a(j);
%     modelt(i,8) = stat(i).phi(j);
%     modelt(i,9) = stat(i).pl(j);
%     if i == 1
%         modelt(i,10) = 0;
%         modelt(i,11) = 0;
%     else
%         modelt(i,10) = stat(i-1).sr(j);
%         modelt(i,11) = stat(i-1).dp(j);
%     end
% end
% [stpl(j).bootmt,stpl(j).rtt,stpl(j).ttt] = forwardHDray2(modelt,param,geom,phaselist,nseg);
% disp([num2str(j) ' of ' num2str(length(stat(1).vp))])
% end
% if ~isempty(modpath)
%     [bootmo,rto,tto] = forwardHDray2(modelo,param,geom,t_ref,phaselist2,nseg2);
% end
fig = figure;%('visible','off')
subplot(1,2,1);
hold on
for i = 1:5
    plot(zz,bootm0(:,i)+i,'k','LineWidth',1.2);
%     for j = 1:length(stpl)
%         plot(time(time<=param(8))-0*param(3),stpl(j).bootmt(time<=param(8),i)+i,'Color',[1 0 0 0.01])
%     end
    plt(i) = plot(zz,bootmr(:,i)+i,'Color',[0 174 0 255/2]./255,'LineWidth',2);
    if ~isempty(modpath)
    plot(zz,bootmo(:,i)+i,'r--','LineWidth',1.1)
    end
end
axis([min(zz) max(zz) 0 6])
xlabel('z in km')
ylabel('Harmonic degrees')
title('modeled')
subplot(1,2,2)
hold on
for i = 6:10
    plot(zz,bootm0(:,i)+i-5,'k','LineWidth',1.2);
%     for j = 1:length(stpl)
%         plot(time(time<=param(8)),stpl(j).bootmt(time<=param(8),i)+i-5,'Color',[1 0 0 0.01])
%     end
    plt2(i-5) = plot(zz,bootmr(:,i)+i-5,'Color',[0 174 0 255/2]./255,'LineWidth',2);
    if ~isempty(modpath)
    plot(zz,bootmo(:,i)+i-5,'r--','LineWidth',1.1)
    end
end
axis([min(zz) max(zz) 0 6])
xlabel('z in km')
title('unmodeled')
print(fig,'-r300','-djpeg','-painters',[savepath '/HarmonicsFinal' addname '.jpg']);
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,'-r0','-dpdf','-vector','-bestfit',[savepath '/HarmonicsFinal' addname '.pdf']);
%print(fig,'-r0','-depsc2','-vector',[savepath '/HarmonicsFinal' addname '.eps']);
if exist('RF','var')
    save([savepath '/statresults' addname '.mat'],'stat','MC','bootm0','RF','bootmr')
else
    save([savepath '/statresults' addname '.mat'],'stat','MC','bootm0','bootmr')
end