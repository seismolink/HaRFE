function MCMCcntrl(opts)

% Bayesian inversion for anisotropic crustal model allowing for dipping
% layers based on a Markov chain Monte Carlo algorithm with Metropolis
% Hastings acceptance criterion
% ---------------------------
% Data adaptive version

% get current folder name
cur_dir = pwd;
[~, df, ~] = fileparts(cur_dir); 

% add all sub directories to search path
addpath(genpath(['../',df]));

if ~isstruct(opts)
    if strcmp(opts,'plot')
        plotallHDMCMCres2;
        return
    end
end

showflag = opts.showflag;

projname = opts.projname; % Name for results and graphics
savepath = [opts.base_dir '/' projname '/'];
if ~exist(savepath,'dir')
    mkdir(savepath);
end

% Reinitialize random number generator
% rng('shuffle');
rng(double(opts.addname))

% Initialize model if opts.init == 1
% ----------------------------------
if opts.init
    load([savepath '/init.mat']);
    for i = 1:length(stat)-1
        if i == 1
            opts.b.hmin(i) =   stat(i).mh-2*stat(i).dh; % thickness
            opts.b.hmax(i) =   stat(i).mh+2*stat(i).dh;
        else
            mh = 0;
            for jj = 1:i
                mh = mh+stat(jj).mh;
            end
            opts.b.hmin(i) =   mh-2*stat(i).dh; % thickness
            opts.b.hmax(i) =   mh+2*stat(i).dh;
        end
        opts.b.vpmin(i) =  stat(i).mvp-2*stat(i).dvp; % p-wave velocity
        opts.b.vpmax(i) =  stat(i).mvp+2*stat(i).dvp;
        % opts.b.vsmin(i) = stat(i).mvp./(stat(i).mvs+2*stat(i).dvs); % vp/vs ratio
        % opts.b.vsmax(i) = stat(i).mvp./(stat(i).mvs-2*stat(i).dvs);
        opts.b.vsmin(i) = stat(i).mvs-2*stat(i).dvs;
        opts.b.vsmax(i) = stat(i).mvs+2*stat(i).dvs;
        opts.b.rhomin(i) = stat(i).mrho-2*stat(i).drho; % density
        opts.b.rhomax(i) = stat(i).mrho+2*stat(i).drho;
        opts.b.amin(i) =   stat(i).ma-2*stat(i).da; % strength of anisotropy in percentage
        opts.b.amax(i) =   stat(i).ma+2*stat(i).da;
        opts.b.phimin(i) = stat(i).mphi-2*stat(i).dphi; % direction of the symmetry axis in degrees
        opts.b.phimax(i) = stat(i).mphi+2*stat(i).dphi;
        opts.b.plmin(i) =  stat(i).mpl-2*stat(i).dpl; % dip of the symmetry axis in degrees
        opts.b.plmax(i) =  stat(i).mpl+2*stat(i).dpl;
        opts.b.srmin(i) = stat(i).msr-2*stat(i).dsr; % strike direction of dipping interface to north
        opts.b.srmax(i) = stat(i).msr+2*stat(i).dsr;
        opts.b.dpmin(i) = 0;%stat(i).mdp-2*stat(i).ddp; % dip of the interface from horizontal
        opts.b.dpmax(i) = stat(i).mdp+2*stat(i).ddp;
        if i == 1
            opts.b.dvp(1) = (stat(length(stat)).mvp-2*stat(length(stat)).dvp)-(stat(length(stat)-1).mvp+2*stat(length(stat)-1).dvp); % velocity contrast between last layer and halfspace
            opts.b.dvp(2) = (stat(length(stat)).mvp+2*stat(length(stat)).dvp)-(stat(length(stat)-1).mvp-2*stat(length(stat)-1).dvp);
        end
    end
    opts.w = ones(opts.nlay-1,10);
    opts.w(2:end,10) = 0;
    
    opts.wi = zeros(size(opts.w));
    opts.wi(:,1) = 0.5;
    opts.wi(:,2) = 0.5;
    opts.wi(:,3) = 0.5;
    opts.wi(:,4) = 0.5;
    opts.wi(:,5) = 0.5;
    opts.wi(:,6) = 0.5;
    opts.wi(:,7) = 0.5;
    opts.wi(:,8) = 0.5;
    opts.wi(:,9) = 0.01;%0.5;
    opts.wi(1,10) = 0.5;
    initflag = 1;
    clear MC RF stat bootm0
else
    initflag = 0;
end

% Initialize parameters for forward calculation
% ------------------------------------
param = opts.param;
% geom = opts.geom;
nlay = opts.nlay;
vpflag = opts.fixvp;
rhflag = opts.fixrh;

disp(num2str(opts.nlay));

time = linspace(0,(param(2)-1).*param(3),param(2))+param(6);

% Load phase list for chosen layered model
% ------------------------------------
if opts.param(1) == 1 
    filename = ['raylist_' num2str(nlay) '_layersPS.txt']; 
else
    filename = ['raylist_' num2str(nlay) '_layers.txt'];
end
[ray] = readphaselist(filename);
segmax = length(ray(end).path);
numph = length(ray);
phaselist = zeros(segmax,2,numph);
for i = 1:length(ray)
    nseg(i) = length(ray(i).path);
    phaselist(1:nseg(i),1,i) = ray(i).path;
    phaselist(1:nseg(i),2,i) = ray(i).type;
end

if opts.param(1) == 1 %&& length(opts.param) < 10
    llflag = 1;
else
    % if opts.param(10) == 1
        llflag = 0;
    % else
    %     llflag = 1;
    % end
end
if llflag
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
else
    nseg2 = nseg;
    phaselist2 = phaselist;
end

% Load data for inversion
% -----------------------
zz = -25:0.25:150;
load(opts.input)
hm(:,zz>param(9)) = [];
errhm(:,zz>param(9)) = [];
zz(zz>param(9)) = [];
hm(:,zz<param(8)) = [];
errhm(:,zz<param(8)) = [];
zz(zz<param(8)) = [];
bootm0 = hm'./max(abs(hm(:)));
ebootm0 = errhm'./max(abs(hm(:)));

if length(opts.param)>9
    rrf = zeros(length(RF),length(time));
    trf = rrf;
    for i = 1:length(RF)
        rrf(i,:) = interp1(RF(i).time,RF(i).rad,time)./max(abs(RF(i).rad));
        trf(i,:) = interp1(RF(i).time,RF(i).tra,time)./max(abs(RF(i).rad));
        geom(i,1) = RF(i).baz;
        geom(i,2) = RF(i).slow;
    end
    rrf(isnan(rrf)) = 0;
    trf(isnan(trf)) = 0;
    rrf(:,time>param(10)) = [];
    trf(:,time>param(10)) = [];
end

% Setup Parameters for the inversion
% ----------------------------------
% minimum and maximum value for each layer, can also be chosen individually
% for each layer above the halfspace
b = opts.b;

% Define free parameters and set fixed parameters of the model
% ------------------------------------------------------------
% 1 = free parameter to search
% 0.[0-99] = fixed parameter to fraction in interval min and max
% -1 = same as previous layer/swap with previous layer (anisotropy)
w = opts.w;

% initial model for free parameters (0 for random initial value)
wi = opts.wi;
n = size(w);
wh = opts.wh;
ws = opts.ws;

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

addname = opts.addname;

% Initialize shape of likelyhood function (Gaussian with width sigt)
roi = 0;
sigt = 0;
%fac = 5/6;%3/5 % Changing this value has a direct influence on the acceptance rate in the Markov chain
fac = 0.33;
Hkamp = 0;
[~,minidx] = min(time.^2);
for i = 1:length(RF)
    [~,idx] = findpeaks(abs(rrf(i,:)));
    idx(idx<minidx) = [];
    ar = rrf(i,idx);
    idx(abs(ar)<0.1.*max(abs(ar))) = [];
    ar(abs(ar)<0.1.*max(abs(ar))) = [];
    [~,idxt] = findpeaks(abs(trf(i,:)));
    idxt(idxt<minidx) = [];
    at = trf(i,idxt);
    idxt(abs(at)<0.1.*max(abs(at))) = [];
    at(abs(at)<0.1.*max(abs(at))) = [];
    RF(i).idx = idx;
    RF(i).idxt = idxt;
    Hkamp = sum(abs(ar))+sum(abs(at))+Hkamp;
end
fnormHk = Hkamp;
sigHk = max(ebootm0(:));
% Hkamp = 1;
% Hkamp2 = 0;%fac./2.*Hkamp;
% while roi<0.5
%     sigt = sigt+0.01;
%     % [L1,fc1] = LeastSquaresLikelyhood_Hkgen(bootm0(:,1:5)',bootm0(:,1:5)',Hkamp,sigt,sigt,wh);
%     % [L2,fc2] = LeastSquaresLikelyhood_Hkgen(bootm0(:,1:5)',bootm0(:,1:5)'.*fac,Hkamp2,sigt,sigt,wh);
%     [L1,fc1] = LeastSquaresLikelyhood_Hkgen(bootm0(zz>0,1:5)',bootm0(zz>0,1:5)',Hkamp,sigt,sigt,wh);
%     [L2,fc2] = LeastSquaresLikelyhood_Hkgen(bootm0(zz>0,1:5)',bootm0(zz>0,1:5)'.*fac,Hkamp2,sigt,sigt,wh);
%     roi = L2./L1;
% end

% no over-writing
A = dir([savepath '/MCMCresult' addname '*.mat']);
nprev = length(A);

% search for free parameters defined in w
N = opts.N; % Number of trials 
Nmc = opts.Nmc; % Number of Markov chains with random starting model (where it is allowed to be random - 0 in wi and 1 in w)
gamma = opts.gamma; % initial step size (decreases adaptive to the acceptance rate)
smoothpar = opts.smoothpar; % parameter for penalizing strong anisotropic models
%maxval = zeros(Nmc,1);
for i = 1:Nmc
    add = [addname '_' num2str(i+nprev)];
    MC.geom = geom;
    MC.param = param;
    MC.nlay = nlay;
    MC.ray = ray;
    MC.w = w;
    MC.wi = wi;
    MC.gamma = gamma;
    MC.N = N;
    MC.b = b;
    MC.initflag = initflag;
    % MC.post = DoMCMCmhc_adap(bootm0(:,1:5),rrf,trf,RF,sigt,b,param,geom,geom2,phaselist2,nseg2,fnormHk,N,gamma,w,wi,wh,showflag,savepath,add,initflag);
    MC.post = DoMCMCmhc_adap(bootm0(:,1:5),rrf,trf,RF,ebootm0,sigHk,b,param,geom,phaselist2,nseg2,fnormHk,N,gamma,smoothpar,w,wi,wh,ws,showflag,savepath,add,vpflag,rhflag,initflag);
    save([savepath '/MCMCresult' add '.mat'],'param','geom','MC')
end

end