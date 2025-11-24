function plotallHDMCMCres

% Bayesian inversion for anisotropic crustal model allowing for dipping
% layers based on a Markov chain Monte Carlo algorithm with Metropolis
% Hastings acceptance criterion
% ---------------------------
% Data adaptive version

NN1 = 5000;
NN2 = 100;

% get current folder name
cur_dir = pwd;
cd ..
programpath = pwd;
[~, df, ~] = fileparts(programpath); 

% add all sub directories to search path
addpath(genpath(['../',df]));
cd(cur_dir)

results_dir = input('Define path to results: ','s');
input_file = input('Define absolute path to input file: ','s');
savepath = results_dir;

addname = 'Final';

A = dir([results_dir '/MCMCresult*.mat']);
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

% Load data for inversion
% -----------------------
load(input_file)
time = linspace(0,(param(2)-1).*param(3),param(2))+param(6);
bootm0 = hm'./max(abs(hm(:)));
time0 = time;
bootm0(end:length(time0),:) = 0;
bootm0(time0>param(8),:) = 0;

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
        if isfield('full',MC(i).post)
        ll = MC(i).post.full.L1';
        res = MC(i).post.full.stat';
        else
        ll = MC(i).post.ll; 
        res = MC(i).post.stat(1:end-1,:); 
        end
    else
        if isfield('full',MC(i).post)
        ll = [MC(i).post.full.L1' ll];
        res = [MC(i).post.full.stat' res];
        else
        ll = [MC(i).post.ll ll]; 
        res = [MC(i).post.stat(1:end-1,:) res];
        end
    end
end
[~,ii] = sort(ll,'descend');
Nres = min(NN1,round(length(ii)./10.*9)); % Number of results accepted from all trials (Nmc*N)
I = ii(1:Nres);
loi = ll(I);

for j = 1:length(ll(I))
    res2(:,:,j) = makemodel(reshape(res(:,I(j)),size(w)),b);
    grvec(:,j) = res(:,I(j));
    n = 1;
    for i = 1:length(res2(:,1,j))
        disp([num2str(j) ' of ' num2str(length(ll(I))) ' - ' num2str(i) ' of ' num2str(length(res2(:,1,j)))]);
        if i == 1
            hr(j,n) = 0;
        else
            hr(j,n) = hr(j,n-1);
        end
        hr(j,n+1) = hr(j,n)+res2(i,1,j)./1000;
        ht(j,n) = res2(i,1,j)./1000;
        ht(j,n+1) = res2(i,1,j)./1000;
        rhor(j,n) = res2(i,2,j)./1000;
        rhor(j,n+1) = res2(i,2,j)./1000;
        vpr(j,n) = res2(i,3,j)./1000;
        vpr(j,n+1) = res2(i,3,j)./1000;
        vsr(j,n) = res2(i,4,j)./1000;
        vsr(j,n+1) = res2(i,4,j)./1000;
        ar(j,n) = res2(i,6,j);
        ar(j,n+1) = res2(i,6,j);
        phir(j,n) = mod(res2(i,8,j),180);
        phir(j,n+1) = mod(res2(i,8,j),180);
        plr(j,n) = res2(i,9,j);
        plr(j,n+1) = res2(i,9,j);
        srr(j,n) = mod(res2(i,10,j),360);
        srr(j,n+1) = mod(res2(i,10,j),360);
        dpr(j,n) = res2(i,11,j);
        dpr(j,n+1) = res2(i,11,j);
        n = n+2;
    end
end
N1 = size(grvec,1);
grveco = grvec;
n = 1;
for i = 1:N1
    if std(grveco(i,:))<eps*1000
        todel(n) = i;
        n = n+1;
        continue
    else
        temp = grvec(i,:);
        temp = temp-min(temp(:));
        temp = temp./max(temp(:));
        grvec(i,:) = temp;
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
Z = linkage(dismat,'complete');
lt = 0; ng = Nres; Ngroups = 3;
while ng(lt==max(lt)) > length(ll(I))/15%min(1000,length(ll)/15)
Ngroups = Ngroups+1;
groups = cluster(Z,'MaxClust',Ngroups);
indx = 1:length(groups);
groupsut = unique(groups);
jj = [];
for i = 1:length(groupsut)
    jj = [jj indx(groupsut(i)==groups)];
    lt(i) = max(loi(groupsut(i)==groups));
    ng(i) = sum(groups==groupsut(i));
end
end

% keyboard
% for i = 1:length(groupsut)
%     ng(i) = sum(groups==groupsut(i));
% end
% [~,ii] = max(ng);
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

alph = ll(I);
alph = alph-min(alph);
alph = alph./max(alph)/length(I)/2;
alph = alph.^0.5;
% alph = alph.^0.2;

% Evaluate results for vp
% -----------------------
fig = figure;
subplot(nlay-1,3,round(1:3:3*(nlay-1)))
hold on
for j = 1:length(I)
    plot(vpr(j,:),hr(j,:),'Color',[0 0 0 alph(j)],'LineWidth',1)
end
set(gca, 'YDir','reverse')
xlabel([parnames{2} ' in [' parunit{2} ']'])
ylabel(['Depth in [' parunit{1} ']'])
for i = 1:nlay-1
    subplot(nlay-1,3,i*3-1)
    hold on
    [Nt,binedges] = histcounts(vpr(:,2*i-1),'BinMethod','auto');
    histogram('BinEdges',binedges,'BinCounts',Nt);
    ymax = max(Nt);
    stat(i).vp = vpr(:,2*i-1);
    stat(i).mvp = vpr(1,2*i-1);%mean(vpr(:,2*i-1));
    stat(i).dvp = 2*sqrt(var(vpr(:,2*i-1)));
    stat(i).nvp = length(vpr(:,2*i-1));
    xstat = linspace(min(binedges),max(binedges),200);
    sigma = stat(i).dvp/2;
    bt = xstat;
    aa = exp(-(bt-stat(i).mvp).^2./(2*sigma*sigma))./(sqrt(2*pi)*sigma).*stat(i).nvp.*(binedges(2)-binedges(1));
    plot(bt,aa,'LineWidth',2)
    ymax = max([ymax max(aa)]);
    %text(stat(i).mvp+(xstat(2)-xstat(1)),max(aa)*1.05,[num2str(round(stat(i).mvp*10)/10) ' ' char(177) ' ' num2str(round(stat(i).dvp*10)/10) ' ' parunit{2}],'VerticalAlignment','bottom')
    plot([stat(i).mvp stat(i).mvp],[0 1.15*ymax],'r-.','LineWidth',1.1)
    set(gca,'YTick',[])
    box('on')
    axis tight
%     if w(i,2) == 1
%     end
    
    subplot(nlay-1,3,i*3)
    hold on
    [Nt,binedges] = histcounts(ht(:,2*i-1),'BinMethod','auto');
    histogram('BinEdges',binedges,'BinCounts',Nt);
    ymax = max(Nt);
    stat(i).h = ht(:,2*i-1);
    stat(i).mh = ht(1,2*i-1);%mean(ht(:,2*i-1));
    stat(i).dh = 2*sqrt(var(ht(:,2*i-1)));
    stat(i).nh = length(ht(:,2*i-1));
    xstat = linspace(min(binedges),max(binedges),200);
    sigma = stat(i).dh/2;
    bt = xstat;
    aa = exp(-(bt-stat(i).mh).^2./(2*sigma*sigma))./(sqrt(2*pi)*sigma).*stat(i).nh.*(binedges(2)-binedges(1));
    plot(bt,aa,'LineWidth',2)
    ymax = max([ymax max(aa)]);
    %text(stat(i).mh+(xstat(2)-xstat(1)),max(aa)*1.05,[num2str(round(stat(i).mh*10)/10) ' ' char(177) ' ' num2str(round(stat(i).dh*10)/10) ' ' parunit{1}],'VerticalAlignment','bottom')
    plot([stat(i).mh stat(i).mh],[0 1.15*ymax],'r-.','LineWidth',1.1)
    set(gca,'YTick',[])
    box('on')
    axis tight
    %     if w(i,1) == 1
%     end
    
end
subplot(nlay-1,3,(nlay-1)*3-1)
xlabel(parnames{2})
subplot(nlay-1,3,(nlay-1)*3)
xlabel(parnames{1})
print(fig,'-r300','-djpeg','-painters',[savepath '/results_vp' addname '.jpg']);

close(fig);


% Evaluate results for vs
% -----------------------
fig = figure;
subplot(nlay-1,3,round(1:3:3*(nlay-1)))
hold on
for j = 1:length(I)
    plot(vsr(j,:),hr(j,:),'Color',[0 0 0 alph(j)],'LineWidth',1)
end
xlabel([parnames{3} ' in [' parunit{3} ']'])
ylabel(['Depth in [' parunit{1} ']'])
set(gca, 'YDir','reverse')
for i = 1:nlay-1
    subplot(nlay-1,3,i*3-1)
    hold on
    [Nt,binedges] = histcounts(vsr(:,2*i-1),'BinMethod','auto');
    histogram('BinEdges',binedges,'BinCounts',Nt);
    ymax = max(Nt);
    stat(i).vs = vsr(:,2*i-1);
    stat(i).mvs = vsr(1,2*i-1);%mean(vsr(:,2*i-1));
    stat(i).dvs = 2*sqrt(var(vsr(:,2*i-1)));
    stat(i).nvs = length(vsr(:,2*i-1));
    xstat = linspace(min(binedges),max(binedges),200);
    sigma = stat(i).dvs/2;
    bt = xstat;
    aa = exp(-(bt-stat(i).mvs).^2./(2*sigma*sigma))./sqrt(pi*2*sigma*sigma)*stat(i).nvs*(binedges(2)-binedges(1));
    plot(bt,aa,'LineWidth',2)
    ymax = max([ymax max(aa)]);
    %text(stat(i).mvs+(xstat(2)-xstat(1)),max(aa)*1.05,[num2str(round(stat(i).mvs*10)/10) ' ' char(177) ' ' num2str(round(stat(i).dvs*10)/10) ' ' parunit{3}],'VerticalAlignment','bottom')
    plot([stat(i).mvs stat(i).mvs],[0 1.15*ymax],'r-.','LineWidth',1.1)
    set(gca,'YTick',[])
    box('on')
    axis tight
%     if w(i,2) == 1
%     end
    
    subplot(nlay-1,3,i*3)
    hold on
    [Nt,binedges] = histcounts(ht(:,2*i-1),'BinMethod','auto');
    histogram('BinEdges',binedges,'BinCounts',Nt);
    ymax = max(Nt);
    xstat = linspace(min(binedges),max(binedges),200);
    sigma = stat(i).dh/2;
    bt = xstat;
    aa = exp(-(bt-stat(i).mh).^2./(2*sigma*sigma))./(sqrt(2*pi)*sigma).*stat(i).nh.*(binedges(2)-binedges(1));
    plot(bt,aa,'LineWidth',2)
    ymax = max([ymax max(aa)]);
    %text(stat(i).mh+(xstat(2)-xstat(1)),max(aa)*1.05,[num2str(round(stat(i).mh*10)/10) ' ' char(177) ' ' num2str(round(stat(i).dh*10)/10) ' ' parunit{1}],'VerticalAlignment','bottom')
    plot([stat(i).mh stat(i).mh],[0 1.15*ymax],'r-.','LineWidth',1.1)
    set(gca,'YTick',[])
    box('on')
    axis tight
    %     if w(i,1) == 1
%     end
    
end
subplot(nlay-1,3,(nlay-1)*3-1)
xlabel(parnames{3})
subplot(nlay-1,3,(nlay-1)*3)
xlabel(parnames{1})
print(fig,'-r300','-djpeg','-painters',[savepath '/results_vs' addname '.jpg']);

close(fig);

% Evaluate results for rho
% -----------------------
fig = figure;
subplot(nlay-1,3,round(1:3:3*(nlay-1)))
hold on
for j = 1:length(I)
    plot(rhor(j,:),hr(j,:),'Color',[0 0 0 alph(j)],'LineWidth',1)
end
xlabel([parnames{4} ' in [' parunit{4} ']'])
ylabel(['Depth in [' parunit{1} ']'])
set(gca, 'YDir','reverse')
for i = 1:nlay-1
    subplot(nlay-1,3,i*3-1)
    hold on
    [Nt,binedges] = histcounts(rhor(:,2*i-1),'BinMethod','auto');
    histogram('BinEdges',binedges,'BinCounts',Nt);
    ymax = max(Nt);
    stat(i).rho = rhor(:,2*i-1);
    stat(i).mrho = rhor(1,2*i-1);%mean(rhor(:,2*i-1));
    stat(i).drho = 2*sqrt(var(rhor(:,2*i-1)));
    stat(i).nrho = length(rhor(:,2*i-1));
    xstat = linspace(min(binedges),max(binedges),200);
    sigma = stat(i).drho/2;
    bt = xstat;
    aa = exp(-(bt-stat(i).mrho).^2./(2*sigma*sigma))./(sqrt(2*pi)*sigma).*stat(i).nrho.*(binedges(2)-binedges(1));
    plot(bt,aa,'LineWidth',2)
    ymax = max([ymax max(aa)]);
    %text(stat(i).mrho+(xstat(2)-xstat(1)),max(aa)*1.05,[num2str(round(stat(i).mrho*10)/10) ' ' char(177) ' ' num2str(round(stat(i).drho*10)/10) ' ' parunit{4}],'VerticalAlignment','bottom')
    plot([stat(i).mrho stat(i).mrho],[0 1.15*ymax],'r-.','LineWidth',1.1)
    set(gca,'YTick',[])
    box('on')
    axis tight
%     if w(i,2) == 1
%     end
    
    subplot(nlay-1,3,i*3)
    hold on
    [Nt,binedges] = histcounts(ht(:,2*i-1),'BinMethod','auto');
    histogram('BinEdges',binedges,'BinCounts',Nt);
    ymax = max(Nt);
    xstat = linspace(min(binedges),max(binedges),200);
    sigma = stat(i).dh/2;
    bt = xstat;
    aa = exp(-(bt-stat(i).mh).^2./(2*sigma*sigma))./(sqrt(2*pi)*sigma).*stat(i).nh.*(binedges(2)-binedges(1));
    plot(bt,aa,'LineWidth',2)
    ymax = max([ymax max(aa)]);
    %text(stat(i).mh+(xstat(2)-xstat(1)),max(aa)*1.05,[num2str(round(stat(i).mh*10)/10) ' ' char(177) ' ' num2str(round(stat(i).dh*10)/10) ' ' parunit{1}],'VerticalAlignment','bottom')
    plot([stat(i).mh stat(i).mh],[0 1.15*ymax],'r-.','LineWidth',1.1)
    set(gca,'YTick',[])
    box('on')
    axis tight
    %     if w(i,1) == 1
%     end
    
end
subplot(nlay-1,3,(nlay-1)*3-1)
xlabel(parnames{4})
subplot(nlay-1,3,(nlay-1)*3)
xlabel(parnames{1})
print(fig,'-r300','-djpeg','-painters',[savepath '/results_rho' addname '.jpg']);

close(fig);


% Evaluate results for a
% -----------------------
fig = figure;
subplot(nlay-1,3,round(1:3:3*(nlay-1)))
hold on
for j = 1:length(I)
    plot(ar(j,:),hr(j,:),'Color',[0 0 0 alph(j)],'LineWidth',1)
end
xlabel([parnames{5} ' in [' parunit{5} ']'])
ylabel(['Depth in [' parunit{1} ']'])
set(gca, 'YDir','reverse')
for i = 1:nlay-1
    subplot(nlay-1,3,i*3-1)
    hold on
    [Nt,binedges] = histcounts(ar(:,2*i-1),'BinMethod','auto');
    histogram('BinEdges',binedges,'BinCounts',Nt);
    ymax = max(Nt);
    stat(i).a = ar(:,2*i-1);
    stat(i).ma = ar(1,2*i-1);%mean(ar(:,2*i-1));
    stat(i).da = 2*sqrt(var(ar(:,2*i-1)));
    stat(i).na = length(ar(:,2*i-1));
    xstat = linspace(min(binedges),max(binedges),200);
    sigma = stat(i).da/2;
    bt = xstat;
    aa = exp(-(bt-stat(i).ma).^2./(2*sigma*sigma))./(sqrt(2*pi)*sigma).*stat(i).na.*(binedges(2)-binedges(1));
    plot(bt,aa,'LineWidth',2)
    ymax = max([ymax max(aa)]);
    %text(stat(i).ma+(xstat(2)-xstat(1)),max(aa)*1.05,[num2str(round(stat(i).ma*10)/10) ' ' char(177) ' ' num2str(round(stat(i).da*10)/10) ' ' parunit{5}],'VerticalAlignment','bottom')
    plot([stat(i).ma stat(i).ma],[0 1.15*ymax],'r-.','LineWidth',1.1)
    set(gca,'YTick',[])
    box('on')
    axis tight
%     if w(i,2) == 1
%     end
    
    subplot(nlay-1,3,i*3)
    hold on
    [Nt,binedges] = histcounts(ht(:,2*i-1),'BinMethod','auto');
    histogram('BinEdges',binedges,'BinCounts',Nt);
    ymax = max(Nt);
    xstat = linspace(min(binedges),max(binedges),200);
    sigma = stat(i).dh/2;
    bt = xstat;
    aa = exp(-(bt-stat(i).mh).^2./(2*sigma*sigma))./(sqrt(2*pi)*sigma).*stat(i).nh.*(binedges(2)-binedges(1));
    plot(bt,aa,'LineWidth',2)
    ymax = max([ymax max(aa)]);
    %text(stat(i).mh+(xstat(2)-xstat(1)),max(aa)*1.05,[num2str(round(stat(i).mh*10)/10) ' ' char(177) ' ' num2str(round(stat(i).dh*10)/10) ' ' parunit{1}],'VerticalAlignment','bottom')
    plot([stat(i).mh stat(i).mh],[0 1.15*ymax],'r-.','LineWidth',1.1)
    set(gca,'YTick',[])
    box('on')
    axis tight
    %     if w(i,1) == 1
%     end
    
end
subplot(nlay-1,3,(nlay-1)*3-1)
xlabel(parnames{5})
subplot(nlay-1,3,(nlay-1)*3)
xlabel(parnames{1})
print(fig,'-r300','-djpeg','-painters',[savepath '/results_a' addname '.jpg']);

close(fig);


% Evaluate results for phi
% -----------------------
fig = figure;
subplot(nlay-1,3,round(1:3:3*(nlay-1)))
hold on
for j = 1:length(I)
    plot(phir(j,:),hr(j,:),'Color',[0 0 0 alph(j)],'LineWidth',1)
end
xlabel([parnames{6} ' in [' parunit{6} ']'])
ylabel(['Depth in [' parunit{1} ']'])
set(gca, 'YDir','reverse')
for i = 1:nlay-1
    subplot(nlay-1,3,i*3-1)
    hold on
    [Nt,binedges] = histcounts(phir(:,2*i-1),'BinMethod','auto');
    histogram('BinEdges',binedges,'BinCounts',Nt);
    ymax = max(Nt);
    
    stat(i).phi = mod(phir(:,2*i-1),180);
    stat(i).nphi = length(phir(:,2*i-1));
    vartemp = mod(phir(:,2*i-1),180);
    vartemp2 = vartemp;
    vartemp2(vartemp2>90) = vartemp2(vartemp2>90)-180;
    dvartemp1 = 2*sqrt(var(vartemp));
    dvartemp2 = 2*sqrt(var(vartemp2));
    
    if dvartemp1<dvartemp2
        stat(i).mphi = vartemp(1);%mean(vartemp);
        stat(i).dphi = dvartemp1;
    else
        stat(i).mphi = mod(vartemp2(1),180);%mod(mean(vartemp2),180);
        stat(i).dphi = dvartemp2;
    end
    
    xstat = linspace(min(binedges),max(binedges),200);
    sigma = stat(i).dphi/2;
    bt = xstat;
    aa = exp(-(bt-stat(i).mphi).^2./(2*sigma*sigma))./(sqrt(2*pi)*sigma).*stat(i).nphi.*(binedges(2)-binedges(1));
    plot(bt,aa,'LineWidth',2)
    ymax = max([ymax max(aa)]);
    %text(stat(i).mphi+(xstat(2)-xstat(1)),max(aa)*1.05,[num2str(round(stat(i).mphi)) ' ' char(177) ' ' num2str(round(stat(i).dphi)) ' ' parunit{6}],'VerticalAlignment','bottom')
    plot([stat(i).mphi stat(i).mphi],[0 1.15*ymax],'r-.','LineWidth',1.1)
    set(gca,'YTick',[])
    box('on')
    axis tight
%     if w(i,2) == 1
%     end
    
    subplot(nlay-1,3,i*3)
    hold on
    [Nt,binedges] = histcounts(ht(:,2*i-1),'BinMethod','auto');
    histogram('BinEdges',binedges,'BinCounts',Nt);
    ymax = max(Nt);
    xstat = linspace(min(binedges),max(binedges),200);
    sigma = stat(i).dh/2;
    bt = xstat;
    aa = exp(-(bt-stat(i).mh).^2./(2*sigma*sigma))./(sqrt(2*pi)*sigma).*stat(i).nh.*(binedges(2)-binedges(1));
    plot(bt,aa,'LineWidth',2)
    ymax = max([ymax max(aa)]);
    %text(stat(i).mh+(xstat(2)-xstat(1)),max(aa)*1.05,[num2str(round(stat(i).mh*10)/10) ' ' char(177) ' ' num2str(round(stat(i).dh*10)/10) ' ' parunit{1}],'VerticalAlignment','bottom')
    plot([stat(i).mh stat(i).mh],[0 1.15*ymax],'r-.','LineWidth',1.1)
    set(gca,'YTick',[])
    box('on')
    axis tight
    %     if w(i,1) == 1
%     end
    
end
subplot(nlay-1,3,(nlay-1)*3-1)
xlabel(parnames{6})
subplot(nlay-1,3,(nlay-1)*3)
xlabel(parnames{1})
print(fig,'-r300','-djpeg','-painters',[savepath '/results_phi' addname '.jpg']);

close(fig);


% Evaluate results for pl
% -----------------------
fig = figure;
subplot(nlay-1,3,round(1:3:3*(nlay-1)))
hold on
for j = 1:length(I)
    plot(plr(j,:),hr(j,:),'Color',[0 0 0 alph(j)],'LineWidth',1)
end
xlabel([parnames{7} ' in [' parunit{7} ']'])
ylabel(['Depth in [' parunit{1} ']'])
set(gca, 'YDir','reverse')
for i = 1:nlay-1
    subplot(nlay-1,3,i*3-1)
    hold on
    [Nt,binedges] = histcounts(plr(:,2*i-1),'BinMethod','auto');
    histogram('BinEdges',binedges,'BinCounts',Nt);
    ymax = max(Nt);
    
    stat(i).pl = plr(:,2*i-1);
    stat(i).mpl = plr(1,2*i-1);%mean(plr(:,2*i-1));
    stat(i).dpl = 2*sqrt(var(plr(:,2*i-1)));
    stat(i).npl = length(plr(:,2*i-1));
    
    xstat = linspace(min(binedges),max(binedges),200);
    sigma = stat(i).dpl/2;
    bt = xstat;
    aa = exp(-(bt-stat(i).mpl).^2./(2*sigma*sigma))./(sqrt(2*pi)*sigma).*stat(i).npl.*(binedges(2)-binedges(1));
    plot(bt,aa,'LineWidth',2)
    ymax = max([ymax max(aa)]);
    %text(stat(i).mpl+(xstat(2)-xstat(1)),max(aa)*1.05,[num2str(round(stat(i).mpl)) ' ' char(177) ' ' num2str(round(stat(i).dpl)) ' ' parunit{7}],'VerticalAlignment','bottom')
    plot([stat(i).mpl stat(i).mpl],[0 1.15*ymax],'r-.','LineWidth',1.1)
    set(gca,'YTick',[])
    box('on')
    axis tight
%     if w(i,2) == 1
%     end
    
    subplot(nlay-1,3,i*3)
    hold on
    [Nt,binedges] = histcounts(ht(:,2*i-1),'BinMethod','auto');
    histogram('BinEdges',binedges,'BinCounts',Nt);
    ymax = max(Nt);
    xstat = linspace(min(binedges),max(binedges),200);
    sigma = stat(i).dh/2;
    bt = xstat;
    aa = exp(-(bt-stat(i).mh).^2./(2*sigma*sigma))./(sqrt(2*pi)*sigma).*stat(i).nh.*(binedges(2)-binedges(1));
    plot(bt,aa,'LineWidth',2)
    ymax = max([ymax max(aa)]);
    %text(stat(i).mh+(xstat(2)-xstat(1)),max(aa)*1.05,[num2str(round(stat(i).mh*10)/10) ' ' char(177) ' ' num2str(round(stat(i).dh*10)/10) ' ' parunit{1}],'VerticalAlignment','bottom')
    plot([stat(i).mh stat(i).mh],[0 1.15*ymax],'r-.','LineWidth',1.1)
    set(gca,'YTick',[])
    box('on')
    axis tight
    %     if w(i,1) == 1
%     end
    
end
subplot(nlay-1,3,(nlay-1)*3-1)
xlabel(parnames{7})
subplot(nlay-1,3,(nlay-1)*3)
xlabel(parnames{1})
print(fig,'-r300','-djpeg','-painters',[savepath '/results_pl' addname '.jpg']);

close(fig);


% Evaluate results for phi
% -----------------------
fig = figure;
subplot(nlay-1,3,round(1:3:3*(nlay-1)))
hold on
for j = 1:length(I)
    plot(srr(j,:),hr(j,:),'Color',[0 0 0 alph(j)],'LineWidth',1)
end
xlabel([parnames{8} ' in [' parunit{8} ']'])
ylabel(['Depth in [' parunit{1} ']'])
set(gca, 'YDir','reverse')
for i = 1:nlay-1
    subplot(nlay-1,3,i*3-1)
    hold on
    [Nt,binedges] = histcounts(srr(:,2*(i+1)-1),'BinMethod','auto');
    histogram('BinEdges',binedges,'BinCounts',Nt);
    ymax = max(Nt);
    
    stat(i).nsr = length(srr(:,2*(i+1)-1));
    stat(i).sr = mod(srr(:,2*(i+1)-1),360);
    vartemp = mod(srr(:,2*(i+1)-1),360);
    vartemp2 = vartemp;
    vartemp2(vartemp2>180) = vartemp2(vartemp2>180)-360;
    dvartemp1 = 2*sqrt(var(vartemp));
    dvartemp2 = 2*sqrt(var(vartemp2));
    
    if dvartemp1<dvartemp2
        stat(i).msr = vartemp(1);%mean(vartemp);
        stat(i).dsr = dvartemp1;
    else
        stat(i).msr = mod(vartemp2(1),360);%mod(mean(vartemp2),360);
        stat(i).dsr = dvartemp2;
    end
    
    xstat = linspace(min(binedges),max(binedges),200);
    sigma = stat(i).dsr/2;
    bt = xstat;
    aa = exp(-(bt-stat(i).msr).^2./(2*sigma*sigma))./(sqrt(2*pi)*sigma).*stat(i).nsr.*(binedges(2)-binedges(1));
    plot(bt,aa,'LineWidth',2)
    ymax = max([ymax max(aa)]);
    %text(stat(i).msr+(xstat(2)-xstat(1)),max(aa)*1.05,[num2str(round(stat(i).msr)) ' ' char(177) ' ' num2str(round(stat(i).dsr)) ' ' parunit{8}],'VerticalAlignment','bottom')
    plot([stat(i).msr stat(i).msr],[0 1.15*ymax],'r-.','LineWidth',1.1)
    set(gca,'YTick',[])
    box('on')
    axis tight
%     if w(i,2) == 1
%     end
    
    subplot(nlay-1,3,i*3)
    hold on
    [Nt,binedges] = histcounts(ht(:,2*i-1),'BinMethod','auto');
    histogram('BinEdges',binedges,'BinCounts',Nt);
    ymax = max(Nt);
    stat(i).mh = ht(1,2*i-1);%mean(ht(:,2*i-1));
    stat(i).dh = 2*sqrt(var(ht(:,2*i-1)));
    stat(i).nh = length(ht(:,2*i-1));
    xstat = linspace(min(binedges),max(binedges),200);
    sigma = stat(i).dh/2;
    bt = xstat;
    aa = exp(-(bt-stat(i).mh).^2./(2*sigma*sigma))./(sqrt(2*pi)*sigma).*stat(i).nh.*(binedges(2)-binedges(1));
    plot(bt,aa,'LineWidth',2)
    ymax = max([ymax max(aa)]);
    %text(stat(i).mh+(xstat(2)-xstat(1)),max(aa)*1.05,[num2str(round(stat(i).mh*10)/10) ' ' char(177) ' ' num2str(round(stat(i).dh*10)/10) ' ' parunit{1}],'VerticalAlignment','bottom')
    plot([stat(i).mh stat(i).mh],[0 1.15*ymax],'r-.','LineWidth',1.1)
    set(gca,'YTick',[])
    box('on')
    axis tight
    %     if w(i,1) == 1
%     end
    
end
subplot(nlay-1,3,(nlay-1)*3-1)
xlabel(parnames{8})
subplot(nlay-1,3,(nlay-1)*3)
xlabel(parnames{1})
print(fig,'-r300','-djpeg','-painters',[savepath '/results_sr' addname '.jpg']);

close(fig);

% Evaluate results for pl
% -----------------------
fig = figure;
subplot(nlay-1,3,round(1:3:3*(nlay-1)))
hold on
for j = 1:length(I)
    plot(dpr(j,:),hr(j,:),'Color',[0 0 0 alph(j)],'LineWidth',1)
end
xlabel([parnames{9} ' in [' parunit{9} ']'])
ylabel(['Depth in [' parunit{1} ']'])
set(gca, 'YDir','reverse')
for i = 1:nlay-1
    subplot(nlay-1,3,i*3-1)
    hold on
    [Nt,binedges] = histcounts(dpr(:,2*(i+1)-1),'BinMethod','auto');
    histogram('BinEdges',binedges,'BinCounts',Nt);
    ymax = max(Nt);
    
    stat(i).dp = dpr(:,2*(i+1)-1);
    stat(i).mdp = dpr(1,2*(i+1)-1);%mean(dpr(:,2*(i+1)-1));
    stat(i).ddp = 2*sqrt(var(dpr(:,2*(i+1)-1)));
    stat(i).ndp = length(dpr(:,2*(i+1)-1));
    
    xstat = linspace(min(binedges),max(binedges),200);
    sigma = stat(i).ddp/2;
    bt = xstat;
    aa = exp(-(bt-stat(i).mdp).^2./(2*sigma*sigma))./(sqrt(2*pi)*sigma).*stat(i).ndp.*(binedges(2)-binedges(1));
    plot(bt,aa,'LineWidth',2)
    ymax = max([ymax max(aa)]);
    %text(stat(i).mdp+(xstat(2)-xstat(1)),max(aa)*1.05,[num2str(round(stat(i).mdp)) ' ' char(177) ' ' num2str(round(stat(i).ddp)) ' ' parunit{9}],'VerticalAlignment','bottom')
    plot([stat(i).mdp stat(i).mdp],[0 1.15*ymax],'r-.','LineWidth',1.1)
    set(gca,'YTick',[])
    box('on')
    axis tight
%     if w(i,2) == 1
%     end
    
    subplot(nlay-1,3,i*3)
    hold on
    [Nt,binedges] = histcounts(ht(:,2*i-1),'BinMethod','auto');
    histogram('BinEdges',binedges,'BinCounts',Nt);
    ymax = max(Nt);
    stat(i).mh = ht(1,2*i-1);%mean(ht(:,2*i-1));
    stat(i).dh = 2*sqrt(var(ht(:,2*i-1)));
    stat(i).nh = length(ht(:,2*i-1));
    xstat = linspace(min(binedges),max(binedges),200);
    sigma = stat(i).dh/2;
    bt = xstat;
    aa = exp(-(bt-stat(i).mh).^2./(2*sigma*sigma))./(sqrt(2*pi)*sigma).*stat(i).nh.*(binedges(2)-binedges(1));
    plot(bt,aa,'LineWidth',2)
    ymax = max([ymax max(aa)]);
    %text(stat(i).mh+(xstat(2)-xstat(1)),max(aa)*1.05,[num2str(round(stat(i).mh*10)/10) ' ' char(177) ' ' num2str(round(stat(i).dh*10)/10) ' ' parunit{1}],'VerticalAlignment','bottom')
    plot([stat(i).mh stat(i).mh],[0 1.15*ymax],'r-.','LineWidth',1.1)
    set(gca,'YTick',[])
    box('on')
    axis tight
    %     if w(i,1) == 1
%     end
    
end
subplot(nlay-1,3,(nlay-1)*3-1)
xlabel(parnames{9})
subplot(nlay-1,3,(nlay-1)*3)
xlabel(parnames{1})
print(fig,'-r300','-djpeg','-painters',[savepath '/results_dp' addname '.jpg']);

close(fig);

stat(nlay).vp = vpr(:,2*nlay-1);
stat(nlay).mvp = vpr(1,2*nlay-1);%mean(vpr(:,2*nlay-1));
stat(nlay).dvp = 2*sqrt(var(vpr(:,2*nlay-1)));
stat(nlay).nvp = length(vpr(:,2*nlay-1));
stat(nlay).h = ht(:,2*nlay-1);
stat(nlay).mh = ht(1,2*nlay-1);%mean(ht(:,2*nlay-1));
stat(nlay).dh = 2*sqrt(var(ht(:,2*nlay-1)));
stat(nlay).nh = length(ht(:,2*nlay-1));
stat(nlay).vs = vsr(:,2*nlay-1);
stat(nlay).mvs = vsr(1,2*nlay-1);%mean(vsr(:,2*nlay-1));
stat(nlay).dvs = 2*sqrt(var(vsr(:,2*nlay-1)));
stat(nlay).nvs = length(vsr(:,2*nlay-1));
stat(nlay).rho = rhor(:,2*nlay-1);
stat(nlay).mrho = rhor(1,2*nlay-1);%mean(rhor(:,2*nlay-1));
stat(nlay).drho = 2*sqrt(var(rhor(:,2*nlay-1)));
stat(nlay).nrho = length(rhor(:,2*nlay-1));
stat(nlay).a = ar(:,2*nlay-1);
stat(nlay).ma = 0;
stat(nlay).da = 0;
stat(nlay).na = length(ar(:,2*nlay-1));
stat(nlay).phi = mod(phir(:,2*nlay-1),180);
stat(nlay).nphi = length(phir(:,2*nlay-1));
stat(nlay).mphi = 0;
stat(nlay).dphi = 0;
stat(nlay).pl = plr(:,2*nlay-1);
stat(nlay).mpl = 0;
stat(nlay).dpl = 0;
stat(nlay).npl = length(plr(:,2*nlay-1));
stat(nlay).nsr = 0;
stat(nlay).sr = 0;
stat(nlay).msr = 0;
stat(nlay).dsr = 0;
stat(nlay).dp = 0;
stat(nlay).mdp = 0;
stat(nlay).ddp = 0;
stat(nlay).ndp = 0;


for i = 1:length(stat)
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

[bootmr,rtr,ttr] = forwardHDray2(modelr,param,geom,phaselist,nseg);
% bootmr(time>param(8),:) = [];


for j = 1:min(NN2,length(stat(1).vp))
    clear modelt
for i = 1:length(stat)
    modelt(i,1) = stat(i).h(j).*1000;
    modelt(i,2) = stat(i).rho(j).*1000;
    modelt(i,3) = stat(i).vp(j).*1000;
    modelt(i,4) = stat(i).vs(j).*1000;
    modelt(i,5) = stat(i).a(j) == 0;
    modelt(i,6) = stat(i).a(j);
    modelt(i,7) = stat(i).a(j);
    modelt(i,8) = stat(i).phi(j);
    modelt(i,9) = stat(i).pl(j);
    if i == 1
        modelt(i,10) = 0;
        modelt(i,11) = 0;
    else
        modelt(i,10) = stat(i-1).sr(j);
        modelt(i,11) = stat(i-1).dp(j);
    end
end
[stpl(j).bootmt,stpl(j).rtt,stpl(j).ttt] = forwardHDray2(modelt,param,geom,phaselist,nseg);
disp([num2str(j) ' of ' num2str(length(stat(1).vp))])
end

fig = figure;%('visible','off')
subplot(1,2,1);
hold on
for i = 1:5
    plot(time(time<=param(8)),bootm0(time<=param(8),i)+i,'b');
    for j = 1:length(stpl)
        plot(time(time<=param(8))-0*param(3),stpl(j).bootmt(time<=param(8),i)+i,'Color',[1 0 0 0.01])
    end
    plt(i) = plot(time(time<=param(8))-0*param(3),bootmr(time<=param(8),i)+i,'r');
end
axis([min(time) param(8) -1 6])
xlabel('t in seconds')
ylabel('Harmonic degrees')
title('modeled')
subplot(1,2,2)
hold on
for i = 6:10
    plot(time(time<=param(8)),bootm0(time<=param(8),i)+i-5,'b');
    for j = 1:length(stpl)
        plot(time(time<=param(8)),stpl(j).bootmt(time<=param(8),i)+i-5,'Color',[1 0 0 0.01])
    end
    plt2(i-5) = plot(time(time<=param(8)),bootmr(time<=param(8),i)+i-5,'r');
end
axis([min(time) param(8) -1 6])
xlabel('t in seconds')
title('unmodeled')
print(fig,'-r300','-djpeg','-painters',[savepath '/HarmonicsFinal' addname '.jpg']);
save([savepath '/statresults' addname '.mat'],'stat','MC','bootm0','bootmr','modelr')

end