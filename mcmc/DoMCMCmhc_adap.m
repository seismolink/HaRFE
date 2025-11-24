function post = DoMCMCmhc_adap(hm,rrf,trf,RF,sigR,sigHk,b,param,geom,phaselist2,nseg2,fnormHk,N,gamma,smoothpar,w,wi,wh,ws,showflag,savepath,add,vpflag,rhflag,initflag)

% sigR(abs(hm)<max(sigR(:))) = max(sigR(:));
% sigR = 2*sigR;
% sigR(sigR<mean(sigR(:))) = mean(sigR(:));
% sigHk = 0.3;
fnorm = max(abs(hm(:,1)));
sigR = sigR./fnorm;
hm = hm./fnorm;
sigR = sqrt(sigR);
sigR = movmean(sigR,10,2);
% sigR(:) = 0.5;
sigHk = 0.5;
wn = w;
wn(w~=1) = 0;
wn2 = w==-1;
% preallocate the sampling point array
n = size(w); % dimensions of the model
mi = length(w(:));
Sample = zeros(mi,N);
% Samplew = zeros(3,1); % the third one is the weight
% time = linspace(0,(param(2)-1).*param(3),param(2))+param(6);
time = param(6):param(3):param(10);
zz = linspace(param(8),param(9),length(hm(:,1)));
% disp(num2str(length(zz)));
% disp(num2str(length(hm(:,1))));
ldepth = [20, 35, 300]; % depths of interfaces
velp = [5.8, 6.5, 8]; % P wave velocity
vels = [3.4, 3.75, 4.5]; % S wave velocity
for i = 1:length(geom(:,1))
    for j = 1:length(zz)
        t_ref{i}(j) = z2t(zz(j),ldepth,velp,vels,geom(i,2));
    end
end

% priors are all uniform distributions > use rand to generate samples
repflag = 1;
while repflag
idx = rand(n);
idx(w~=1) = w(w~=1);
idx(wi~=0) = wi(wi~=0);
m = makemodel(idx,b);
dh = m(1:end-1,1);
if min(dh) >= 1000
    repflag = 0;
end
end
if ~initflag
fm = max(max(abs(b.amax)),max(abs(b.amin)))./((max(max(abs(b.amax)),max(abs(b.amin)))-10).*(m(:,1)./1000).^-0.25+10);
fm(fm<1) = (m(fm<1,1)./1000).^-2;
m(:,6:7) = m(:,6:7)./fm;
end
Samplew(1:mi,1) = idx(:); Samplew(mi+1,1) = 1;

allsample = zeros(N,length(idx(:)));
allL = zeros(N,1);
allL2 = zeros(N,1);
allLo = zeros(N,1);
allalph = zeros(N,1);
allalph2 = zeros(N,1);
alllim = zeros(N,1);
allHk = zeros(N,1);
allrms = zeros(N,1);

% total model prior density function
% calculate first forward model
[hdc,~,~,wr,wt] = forwardHDray2(m,param,geom,t_ref,phaselist2,nseg2);
hdc = hdc./max(abs(hdc(:,1)));
M1 = size(rrf);
M2 = size(wr);
Hkamp = 0;
for i = 1:M1(1)
    %nn = M2(1)-M1(1)+i;
    nn = i;
    Hkamp = sum(rrf(i,RF(i).idx).*wr(nn,RF(i).idx))+sum(trf(i,RF(i).idxt).*wt(nn,RF(i).idxt))+Hkamp;
end
Hkamp = Hkamp./fnormHk;
if Hkamp > 1
    disp('Hkamp > 1')
    Hkamp(Hkamp>1) = 1;
end
if Hkamp < 0
    Hkamp = 0;
end

% likelyhood function
% [Lold,~] = LeastSquaresLikelyhood_Hkgen(hm',hdc(:,1:5)',Hkamp,sigR,sigT,wh);
% [Lold,~] = LeastSquaresLikelyhood_Hkgen_lowA(hm',hdc(:,1:5)',Hkamp,sigR,sigT,smoothpar,m(1:end-1,6),m(1:end-1,1),wh);
[Lold,~,logLold] = LeastSquaresLikelyhood_Hkgen_lowA(hm(zz>0,1:5)',hdc(zz>0,1:5)',Hkamp,sigR(zz>0,1:5)',sigHk,wh);
ll(1) = logLold;

N0 = 0;
ni = 1;

if showflag
fig = figure;%('visible','off')
subplot(1,3,1);
hold on
hdc0 = hdc;
for i = 1:5
    plot(zz,hm(:,i)+i,'b','LineWidth',1.2);
    plt0(i) = plot(zz,hdc0(:,i)+i,'k');
    plt(i) = plot(zz,hdc(:,i)+i,'r');
end
axis([min(zz) max(zz) 0 6])
xlabel('z in km')
ylabel('Harmonic degrees')
title('modeled')
subplot(1,3,2)
hold on
M = size(rrf);
wr0 = wr;
for i = 1:M(1)
    plot(time,rrf(i,:)./max(abs(rrf(i,:)))+i,'b')
    plot([time(RF(i).idx)' time(RF(i).idx)']',[-1 1]'*ones(1,length(RF(i).idx))+i,'g')
    plt20(i) = plot(time,wr0(i,:)+i,'k');
    plt2(i) = plot(time,wr(i,:)+i,'r');
end
axis([min(time) param(10) 0 M(1)+1])
xlabel('t in seconds')
title('Radial')
subplot(1,3,3)
hold on
nn = 1;
for i = 1:length(m(:,1))
    if i == 1
        hr(nn) = 0;
    else
        hr(nn) = hr(nn-1);
    end
    hr(nn+1) = hr(nn)+m(i,1);
    rhor(nn) = m(i,2);
    rhor(nn+1) = m(i,2);
    vpr(nn) = m(i,3);
    vpr(nn+1) = m(i,3);
    vsr(nn) = m(i,4);
    vsr(nn+1) = m(i,4);
    ar(nn) = m(i,6);
    ar(nn+1) = m(i,6);
    phir(nn) = m(i,7);
    phir(nn+1) = m(i,7);
    plr(nn) = m(i,8);
    plr(nn+1) = m(i,8);
    srr(nn) = m(i,9);
    srr(nn+1) = m(i,9);
    dpr(nn) = m(i,10);
    dpr(nn+1) = m(i,10);
    nn = nn+2;
end
vpr0 = vpr;
vsr0 = vsr;
hr0 = hr;
ar0 = ar;
plt0a = plot(vpr0,hr0,'k');
plt0b = plot((1+abs(ar0)./100).*vpr0,hr0,'k-.');
plt0c = plot((1-abs(ar0)./100).*vpr0,hr0,'k-.');
plt0d = plot(vsr0,hr0,'k');
plt0e = plot((1+abs(ar0)./100).*vsr0,hr0,'k-.');
plt0f = plot((1-abs(ar0)./100).*vsr0,hr0,'k-.');
plt3a = plot(vpr,hr,'r');
plt3b = plot((1+abs(ar)./100).*vpr,hr,'r-.');
plt3c = plot((1-abs(ar)./100).*vpr,hr,'r-.');
plt3d = plot(vsr,hr,'b');
plt3e = plot((1+abs(ar)./100).*vsr,hr,'b-.');
plt3f = plot((1-abs(ar)./100).*vsr,hr,'b-.');
drawnow;
end

sampleCount = 2;
% Metropolis-Hastings loop
while sampleCount <= N
    tic
    % Generate a new point from the current point from the proposal
    % distribution 
    %%% replace with new proposal draw
    % p = idx+gamma*randn(n).*wn;
    % if ~initflag
    % p(:,6) = mod(p(:,6),1);
    % p(:,8) = mod(p(:,8),1);
    % end
    % p(p>1) = 1-mod(p(p>1),1);
    % p(p<0) = 1-mod(p(p<0),1);
    % p(wn2) = -1;
    % while sum(abs(p(~wn2))>1) > 0 || sum(p(~wn2)<0) || any(p(ws==-1,5)>0.5) || any(p(ws==1,5)<0.5)
    %     p = idx+gamma*randn(n).*wn;
    %     if ~initflag
    %     p(:,6) = mod(p(:,6),1);
    %     p(:,8) = mod(p(:,8),1);
    %     end
    %     p(p>1) = 1;
    %     p(p<0) = 0;
    %     p(wn2) = -1;
    % end
    % if ~vpflag
    %     if sum(wn(:,2)) == 0
    %         dbb = (p(:,3)-idx(:,3));%./p(:,3);
    %         p(:,2) = p(:,2)+dbb.*0.8;%.*p(:,2);
    %         %disp(num2str(p(:,2)))
    %     end
    % end
    % if ~rhflag
    %     if sum(wn(:,4)) == 0
    %         dbb = (p(:,3)-idx(:,3));%./p(:,3);
    %         p(:,4) = p(:,4)+dbb.*0.4;%.*p(:,4);
    %         %disp(num2str(p(:,4)))
    %     end
    % end
    % % if sampleCount<N/8
    % %     p(:,5:end) = idx(:,5:end);
    % % end
    % % Compute the target density
    % pm = makemodel(p,b);
    %%% replace with the following while loop
    repflag = 1;
    while repflag
        p = idx+gamma*randn(n).*wn;
        if ~initflag
            p(:,6) = mod(p(:,6),1);
            p(:,8) = mod(p(:,8),1);
        end
        p(p>1) = 1-mod(p(p>1),1);
        p(p<0) = 1-mod(p(p<0),1);
        p(wn2) = -1;

        if ~vpflag
            if sum(wn(:,2)) == 0
                dbb = (p(:,3)-idx(:,3));%./p(:,3);
                p(:,2) = p(:,2)+dbb.*0.8;%.*p(:,2);
                %disp(num2str(p(:,2)))
            end
        end
        if ~rhflag
            if sum(wn(:,4)) == 0
                dbb = (p(:,3)-idx(:,3));%./p(:,3);
                p(:,4) = p(:,4)+dbb.*0.4;%.*p(:,4);
                %disp(num2str(p(:,4)))
            end
        end

        pm = makemodel(p,b);
        dh = pm(1:end-1,1);
        if min(dh) < 1000 || sum(abs(p(~wn2))>1) > 0 || sum(p(~wn2)<0) || any(p(ws==-1,5)>0.5) || any(p(ws==1,5)<0.5)
            repflag = 1;
        else
            repflag = 0;
        end
    end
    % disp([num2str(pm(2,7)) ' ' num2str(pm(2,8)) ' ' num2str(pm(2,9))])
    if ~initflag
    fm = max(max(abs(b.amax)),max(abs(b.amin)))./((max(max(abs(b.amax)),max(abs(b.amin)))-10).*(pm(:,1)./1000).^-0.25+10);
    fm(fm<1) = (pm(fm<1,1)./1000).^-2;
    pm(:,6:7) = pm(:,6:7)./fm;
    end
    try
        [hdc,~,~,wr,wt] = forwardHDray2(pm,param,geom,t_ref,phaselist2,nseg2);
        hdc = hdc./max(abs(hdc(:,1)));
        if sum(hdc(:)) == 0
            continue
        end
        Hkamp = 0;
        %n = 0;
        for i = 1:M1(1)
            %nn = M2(1)-M1(1)+i;
            nn = i;
            Hkamp = sum(rrf(i,RF(i).idx).*wr(nn,RF(i).idx))+sum(trf(i,RF(i).idxt).*wt(nn,RF(i).idxt))+Hkamp;
        end
        Hkamp = Hkamp./fnormHk;
        if Hkamp > 1
            disp('Hkamp > 1')
            Hkamp(Hkamp>1) = 1;
        end
        if Hkamp < 0
            Hkamp = 0;
        end
    catch
        continue
    end
    % [Lnew,~] = LeastSquaresLikelyhood_Hkgen(hm',hdc(:,1:5)',Hkamp,sigR,sigT,wh);
    % [Lnew2,~] = LeastSquaresLikelyhood_Hkgen(hm',hdc(:,1:5)',0,sigR,sigT,wh);
    % [Lnew,~] = LeastSquaresLikelyhood_Hkgen_lowA(hm',hdc(:,1:5)',Hkamp,sigR,sigT,smoothpar,pm(1:end-1,6),pm(1:end-1,1),wh);
    % [Lnew2,~] = LeastSquaresLikelyhood_Hkgen_lowA(hm',hdc(:,1:5)',0,sigR,sigT,smoothpar,pm(1:end-1,6),pm(1:end-1,1),wh);
    [Lnew,~,logLnew] = LeastSquaresLikelyhood_Hkgen_lowA(hm(zz>0,1:5)',hdc(zz>0,1:5)',Hkamp,sigR(zz>0,1:5)',sigHk,wh);
    [Lnew2,~,logLnew2] = LeastSquaresLikelyhood_Hkgen_lowA(hm(zz>0,1:5)',hdc(zz>0,1:5)',0,sigR(zz>0,1:5)',sigHk,wh);
    % Compute the proposal
    % qmp = proposaldensity(idx,p,gamma); qpm = proposaldensity(p,idx,gamma);
    [qmp,logqmp] = proposaldensity(wi,p,smoothpar); [qpm,logqpm] = proposaldensity(wi,idx,smoothpar);
    % Compute the acceptance ratio, using log to avoid underflow and
    % overflow
    % alpha = log(Lnew*qpm)-log(Lold*qmp);
    % alpha2 = log(Lnew2*qpm)-log(Lold*qmp);
    alpha = (logLnew+logqmp)-(logLold+logqpm);
    alpha2 = (logLnew2+logqmp)-(logLold+logqpm);
    % alpha = (logLnew)-(logLold);
    % alpha2 = (logLnew2)-(logLold);

    if showflag
    subplot(1,3,1)
    delete(plt)
    for i = 1:5
        plt(i) = plot(zz,hdc(:,i)+i,'r');
    end
    drawnow
    subplot(1,3,2)
    delete(plt2)
    for i = 1:M(1)
        nn = M2(1)-M1(1)+i;
        plt2(i) = plot(time,wr(nn,:)+i,'r');
    end
    drawnow
    subplot(1,3,3)
    hold on
    nn = 1;
    for i = 1:length(pm(:,1))
        if i == 1
            hr(nn) = 0;
        else
            hr(nn) = hr(nn-1);
        end
        hr(nn+1) = hr(nn)+pm(i,1);
        rhor(nn) = pm(i,2);
        rhor(nn+1) = pm(i,2);
        vpr(nn) = pm(i,3);
        vpr(nn+1) = pm(i,3);
        vsr(nn) = pm(i,4);
        vsr(nn+1) = pm(i,4);
        ar(nn) = pm(i,6);
        ar(nn+1) = pm(i,6);
        phir(nn) = pm(i,7);
        phir(nn+1) = pm(i,7);
        plr(nn) = pm(i,8);
        plr(nn+1) = pm(i,8);
        srr(nn) = pm(i,9);
        srr(nn+1) = pm(i,9);
        dpr(nn) = pm(i,10);
        dpr(nn+1) = pm(i,10);
        nn = nn+2;
    end
    delete(plt3a); delete(plt3b); delete(plt3c); delete(plt3d); delete(plt3e); delete(plt3f);
    plt3a = plot(vpr,hr,'r');
    plt3b = plot((1+abs(ar)./100).*vpr,hr,'r-.');
    plt3c = plot((1-abs(ar)./100).*vpr,hr,'r-.');
    plt3d = plot(vsr,hr,'b');
    plt3e = plot((1+abs(ar)./100).*vsr,hr,'b-.');
    plt3f = plot((1-abs(ar)./100).*vsr,hr,'b-.');
    drawnow;
    end
    
    % if initflag
    %     aplim = 0.01;
    % else
    %     aplim = 0.08;
    % end
    % rx = (rand+(1-aplim)./(1-(1-aplim)))./(1+(1-aplim)./(1-(1-aplim)));
    % alpchck = log(rx);
    alpchck = log(rand);
    while alpchck<-0.05
        alpchck = log(rand);
    end
    % temp1 = hm(:,1:5);
    % temp2 = hdc(:,1:5);
    temp1 = hm(zz>0,1:5);
    temp2 = hdc(zz>0,1:5);
    disp(['Lold = ' num2str(logLold) ' Lnew = ' num2str(logLnew) ' Lnew2 = ' num2str(logLnew2)])
    disp(['Hk = ' num2str(Hkamp) ' rms = ' num2str(sqrt(mean((temp1(:)-temp2(:)).^2))) ' alpha = ' num2str(alpha) ' alpha2 = ' num2str(alpha2) ' lim = ' num2str(alpchck)]);
    allsample(sampleCount,:) = p(:);
    allL(sampleCount) = logLnew;
    allL2(sampleCount) = logLnew2;
    allLo(sampleCount) = logLold;
    allalph(sampleCount) = alpha;
    allalph2(sampleCount) = alpha2;
    alllim(sampleCount) = alpchck;
    allHk(sampleCount) = Hkamp;
    allrms(sampleCount) = sqrt(mean((temp1(:)-temp2(:)).^2));
    if (alpha > alpchck)
        % accept the new sampling point
        idx = p;
        Lold = Lnew;
        logLold = logLnew;
        Samplew = [Samplew,ones(mi+1,1)];
        Samplew(1:mi,end) = idx(:);
        ll = [ll logLold];

        if showflag
        % replot reference model
        subplot(1,3,1)
        hdc0 = hdc;
        delete(plt0)
        for i = 1:5
            plt0(i) = plot(zz,hdc0(:,i)+i,'k');
        end
        drawnow
        subplot(1,3,2)
        wr0 = wr;
        delete(plt20)
        for i = 1:M(1)
            nn = M2(1)-M1(1)+i;
            plt20(i) = plot(time,wr0(nn,:)+i,'k');
        end
        vpr0 = vpr;
        vsr0 = vsr;
        hr0 = hr;
        ar0 = ar;
        subplot(1,3,3)
        hold on
        delete(plt0a); delete(plt0b); delete(plt0c); delete(plt0d); delete(plt0e); delete(plt0f);
        plt0a = plot(vpr0,hr0,'k');
        plt0b = plot((1+abs(ar0)./100).*vpr0,hr0,'k-.');
        plt0c = plot((1-abs(ar0)./100).*vpr0,hr0,'k-.');
        plt0d = plot(vsr0,hr0,'k');
        plt0e = plot((1+abs(ar0)./100).*vsr0,hr0,'k-.');
        plt0f = plot((1-abs(ar0)./100).*vsr0,hr0,'k-.');
        end
    else
        % increase the weight if stay put
        Samplew(mi+1,end) = Samplew(mi+1,end) + 1;
    end
    disp([num2str(sampleCount) ' steptime: ' num2str(toc) ' s'])
    Sample(:,sampleCount) = idx(:);
    sampleCount = sampleCount+1;
    if mod(sampleCount,round(N/10))==0
        if sampleCount>N/4
            aratetemp = (size(Samplew,2)-N0)/N;
            if aratetemp < 0.2
                gamma = gamma/5*3;
                if gamma <= 0.02
                    gamma = 0.02;
                    disp('Step size decreases no more.')
                end
                [~,ii] = max(ll);
                idx = reshape(Samplew(1:mi,ii),size(idx)) ;
            end
        end
        N0 = size(Samplew,2);
    end
end

% average acceptance rate
arate = size(Samplew,2)/N;

N2 = size(Samplew,2);
indx = 1:N2;
i0 = indx(ll>(max(ll)/10));
for i = 1:length(i0)%round(N2/4):N2
    mm = makemodel(reshape(Samplew(1:mi,i0(i)),size(idx)),b);
    if ~initflag
    fm = max(max(abs(b.amax)),max(abs(b.amin)))./((max(max(abs(b.amax)),max(abs(b.amin)))-10).*(mm(:,1)./1000).^-0.25+10);
    fm(fm<1) = (mm(fm<1,1)./1000).^-2;
    mm(:,6:7) = mm(:,6:7)./fm;
    end
    post.para(:,:,ni) = mm;
    post.L(ni) = ll(i0(i));
    ni = ni+1;
end
post.stat = Samplew;
post.ll = ll;
post.arate = arate;
post.ii = i0;
post.full.stat = allsample;
post.full.L1 = allL;
post.full.L2 = allL2;
post.full.Lo = allLo;
post.full.alpha1 = allalph;
post.full.alpha2 = allalph2;
post.full.lim = alllim;
post.full.Hk = allHk;
post.full.rms = allrms;
disp(['Number of accepted solutions: ' num2str(ni-1)])
if showflag
print(fig,[savepath '\FitRFMCMC_' add '.jpg'],'-djpeg','-r600')
close(fig);
end
end