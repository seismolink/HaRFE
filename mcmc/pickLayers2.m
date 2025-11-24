clear all
close all

filedir = input('Define path to results: ','s');
statlist = input('Define absolute path to station file: ','s');
nlay = input('Define number of layers for initial model guess: ');

% load station information
fid = fopen(statlist,'rt');
D = textscan(fid,'%f %f %f %f %s %s','Delimiter',',','Headerlines',1);
fclose(fid);
% station = C{1,1};
lat = D{1,2};
lon = D{1,3};
ele = D{1,4};
statall = D{1,5};
net = D{1,6};
% indx = 1:length(stat);

A = dir([filedir '/*.mat']);
% Define starting point of the profile
disp('The following stations are available as reference: ')
m = 1;
indx = 1:length(statall);
for i = 1:length(A)
    bla = indx(strcmp(A(i).name(1:end-7),statall));
    if sum(bla) > 0
        latf(m) = lat(bla);
        lonf(m) = lon(bla);
        elef(m) = ele(bla);
        statf{m} = statall{bla};
        netf{m} = net{bla};
        disp([num2str(m) ' ' statf{m} ' ' netf{m} ' ' num2str(latf(m)) ' ' num2str(lonf(m))])
        m = m+1;
    end
end
ii = input('Station identifier for reference point: ');
disp(['Station ' statf{ii} ' selected as reference.']);
latref = latf(ii);
lonref = lonf(ii);

stat = statf;
lat = latf;
lon = lonf;
ele = elef;
net = netf;
dd = sqrt((latref-lat).^2+(lonref-lon).^2);
[dd,I] = sort(dd);
stat = stat(I);
lat = lat(I);
lon = lon(I);
ele = ele(I);
net = net(I);

% load harmonics from path and station list
indx = 1:length(stat);
for ii = 1:length(stat)
    A = dir([filedir '/' stat{ii} '*.mat']);
    % prepare folders for data
%     ind = strfind(A(ii).name,'.');
%     stat{ii} = A(ii).name(1:ind(1)-2);
    
    load([filedir '/' A(1).name]);
    
    hm(6:end,:) = [];
    bootm0 = hm'./max(abs(hm(:)));
    hdt = bootm0(:);
    
%     fun(ii,:) = bootm0(:);
    amp1(ii,:) = bootm0(:,1)./max(abs(hdt(:)));
    amp2(ii,:) = bootm0(:,2)./max(abs(hdt(:)));
    amp3(ii,:) = bootm0(:,3)./max(abs(hdt(:)));
    amp4(ii,:) = bootm0(:,4)./max(abs(hdt(:)));
    amp5(ii,:) = bootm0(:,5)./max(abs(hdt(:)));
    fun(ii,:) = [amp1(ii,:)/max(abs(amp1(ii,:))) amp2(ii,:)/max(abs(amp2(ii,:))) amp3(ii,:)/max(abs(amp3(ii,:))) amp4(ii,:)/max(abs(amp4(ii,:))) amp5(ii,:)/max(abs(amp5(ii,:))) dd(ii)./max(dd)];
end


jjn = 1:length(latf);

% preparation of parameters
amp1 = amp1(jjn,:);
amp2 = amp2(jjn,:);
amp3 = amp3(jjn,:);
amp4 = amp4(jjn,:);
amp5 = amp5(jjn,:);
% ampa = (abs(amp1)./max(abs(amp1(:))))+(abs(amp2)./max(abs(amp2(:))))+(abs(amp3)./max(abs(amp3(:))))+(abs(amp4)./max(abs(amp4(:))))+(abs(amp5)./max(abs(amp5(:))));
ampa = 2.*(abs(amp1)./max(abs(amp1(:))))...
    +(abs(amp2)./max(max(abs(amp2(:))),max(abs(amp3(:)))))...
    +(abs(amp3)./max(max(abs(amp3(:))),max(abs(amp3(:)))))...
    +(abs(amp4)./max(max(abs(amp4(:))),max(abs(amp5(:)))))...
    +(abs(amp5)./max(max(abs(amp5(:))),max(abs(amp5(:)))));
ampa = (ampa./max(abs(ampa(:)))).^0.1;
% ampa = (abs(amp1)./max(abs(amp1(:)))).^0.25+(abs(amp2)./max(abs(amp2(:)))).^0.25+(abs(amp3)./max(abs(amp3(:)))).^0.25+(abs(amp4)./max(abs(amp4(:)))).^0.25+(abs(amp5)./max(abs(amp5(:)))).^0.25;
% ampa = smoothdata(ampa,2);
ampa = ampa./max(abs(ampa),[],2);
ampa2 = (abs(amp1)./max(abs(amp1(:)))).^0.25+(abs(amp2)./max(abs(amp2(:)))).^0.25+(abs(amp3)./max(abs(amp3(:)))).^0.25;
ampa2 = smoothdata(ampa2,2);
ampa2 = ampa2./max(abs(ampa2),[],2);
lat = lat(jjn);
lon = lon(jjn);
dd = dd(jjn);
stat = stat(jjn);
[~,iraw] = sort(dd);
w = [1 1 1 1 1];
xvec = indx;
zvec = -25:0.25:150;
%zvec = -2:0.025:12;
xlim = [xvec(1)-0.5 xvec(end)+0.5];
zlim = [zvec(1) zvec(end)];
selmat = zeros(size(amp1));
mohovec = zeros(length(indx),1);
M = size(selmat);

sm = "hexagram";
cm = [1 1 1];
ss = ["o","square","diamond","^","v",">","<","pentagram","+","*",".","x","_","|"];
cs = [   0         0    1.0000;
         0    0.7008         0;
    1.0000         0    0.4444;
    0.5873         0         0;
    0.5873    1.0000         0;
    0.1587    0.5512    1.0000;
         0         0    0.4762;
    1.0000    0.6772    0.0476;
    0.8413         0    1.0000;
    0.5873    0.4882    0.4286;
    1.0000    0.6142    0.9683;
    0.1111    0.8583    0.8730;
    0.1746    0.0866         0;
    0.4603    0.0709    0.6190;
         0    0.3701         0];
   
 % create red waveform
nk = round(256/12);
nw = round(nk/3);
r = zeros(256, 1);
r(1:nk) = linspace(0.6,1,length(1:nk));
r(nk+1:128-nw) = 1;
r(128-nw+1:128+nw) = 1;
r(129+nw:256-nk) = linspace(1,0,length(129+nw:256-nk));
r(256-(nk-1):256) = 0;
% Create green waveform:
g = zeros(256, 1);
g(1:nk) = 0;
g(nk+1:128-nw) = linspace(0,1,length(nk+1:128-nw));
g(128-nw+1:128+nw) = 1;
g(129+nw:256-nk) = linspace(1,0,length(129+nw:256-nk));
g(256-(nk-1):256) = 0;
% Create blue waveform:
b = zeros(256, 1);
b(1:nk) = 0;
b(nk+1:128-nw) = linspace(0,1,length(nk+1:128-nw));
b(129-nw:128+nw) = 1;
b(129-nw:256-nk) = 1;
b(256-(nk-1):256) = linspace(1,0.6,length(256-(nk-1):256));
% % Now that the individual color channel mappings have been created,
% % stitch the red, green, and blue vectors together to create the 256-by-3 colormap matrix.
customColorMap = [r, g, b];

% produce initial guess
ix = 1:length(zvec);
temp = ampa(iraw,:);
tempm = amp1(iraw,:);
ampa3 = zeros(size(ampa));
amp1b = zeros(size(ampa));
for i = 1:length(iraw)
    idx = max(1,i-2):min(length(iraw),i+2);
    temp2 = sum(temp(idx,:),1)+temp(i,:);
    temp2m = sum(tempm(idx,:),1);%+tempm(i,:);
    indxminm = ix(zvec>=2.5)+1;
    indxminm = indxminm(1);
    [pks,locs] = findpeaks(temp2m);
    pks(locs<=indxminm) = [];
    locs(locs<=indxminm) = [];
    [~,I] = sort(pks,'descend');
    pks = pks(I); locs = locs(I);
    locsoi = locs(1);
    xn = ones(size(locsoi)).*i;
    iz = locsoi;
    mohovec(iraw(xn)) = zvec(iz);
    
    [pks,locs] = findpeaks(temp2);
    indxmin = ix(zvec==0)+1;
    pks(locs<=indxmin) = [];
    locs(locs<=indxmin) = [];
    pks(locs>=iz-round(0.5./0.025)) = [];
    locs(locs>=iz-round(0.5./0.025)) = [];
    [~,I] = sort(pks,'descend');
    pks = pks(I); locs = locs(I);
    locsoi = locs(1:min(length(locs),nlay-1));
    xn = ones(size(locsoi)).*i;
    iz = locsoi;
    selmat(iraw(i),iz) = 1;
    ampa3(iraw(i),:) = temp2.^2;
    amp1b(iraw(i),:) = temp2m;
end

finished_flag = 0;
while ~finished_flag
    anyinput = input('PI>>','s');
    if strcmp(anyinput,'quit')
        finished_flag = 1;
        continue
    end

    if strcmp(anyinput,'plotraw') || strcmp(anyinput,'pr')
        if ~exist('figmain','var')
            figmain = figure('units','normalized','outerposition',[0 0.05 1 0.95]);
        else
            figure(figmain);
        end
        hold on
        imagesc(xvec,zvec,ampa(iraw,:)')
        colormap(jet);
        caxis([min(ampa(:)) max(ampa(:))])
        set(gca,'Ydir','reverse')
        plot([xvec(1)-0.5 xvec(end)+0.5],[0 0],'k','LineWidth',1.1)
        axis([xlim zlim])
%         plot([xvec(ig)' xvec(ig)']'+0.5,[zvec(1).*ones(size(xvec(ig)')) zvec(end).*ones(size(xvec(ig)'))]','r','LineWidth',1.1)
%         plot(xvec(indxn),zeros(length(indxn),1),'w*','MarkerSize',8,'LineWidth',1.1)
        n = 1;
        for i = 1:length(indx)
            temp = selmat(iraw(i),:);
            it = find(temp);
            if ~isempty(it)
                for j = 1:length(it)
                    pl(n) = plot(i,zvec(it(j)),ss(j),'Color','w','MarkerFaceColor',cs(j,:),'MarkerSize',6);
                    n = n+1;
                end
            end
        end
        mpl = find(mohovec);
        if ~isempty(mpl)
            plm = plot(mpl,mohovec(iraw(mpl)),sm,'Color','w','MarkerFaceColor',cm,'MarkerSize',6);
        end
        hold off
        opt = -1;
        continue
    end

    if strcmp(anyinput,'plotrawb') || strcmp(anyinput,'prb')
        if ~exist('figmain','var')
            figmain = figure('units','normalized','outerposition',[0 0.05 1 0.95]);
        else
            figure(figmain);
        end
        hold on
        imagesc(xvec,zvec,ampa2(iraw,:)')
        colormap(jet);
        caxis([min(ampa2(:)) max(ampa2(:))])
        set(gca,'Ydir','reverse')
        plot([xvec(1)-0.5 xvec(end)+0.5],[0 0],'k','LineWidth',1.1)
        axis([xlim zlim])
%         plot([xvec(ig)' xvec(ig)']'+0.5,[zvec(1).*ones(size(xvec(ig)')) zvec(end).*ones(size(xvec(ig)'))]','r','LineWidth',1.1)
%         plot(xvec(indxn),zeros(length(indxn),1),'w*','MarkerSize',8,'LineWidth',1.1)
        n = 1;
        for i = 1:length(indx)
            temp = selmat(iraw(i),:);
            it = find(temp);
            if ~isempty(it)
                for j = 1:length(it)
                    pl(n) = plot(i,zvec(it(j)),ss(j),'Color','w','MarkerFaceColor',cs(j,:),'MarkerSize',6);
                    n = n+1;
                end
            end
        end
        mpl = find(mohovec);
        if ~isempty(mpl)
            plm = plot(mpl,mohovec(iraw(mpl)),sm,'Color','w','MarkerFaceColor',cm,'MarkerSize',6);
        end
        hold off
        opt = -1;
        continue
    end

    if strcmp(anyinput,'plotall') || strcmp(anyinput,'pa')
        if ~exist('figmain','var')
            figmain = figure('units','normalized','outerposition',[0 0.05 1 0.95]);
        else
            figure(figmain);
        end
        hold on
        imagesc(xvec,zvec,ampa')
        colormap(jet);
        caxis([min(ampa(:)) max(ampa(:))])
        set(gca,'Ydir','reverse')
        plot([xvec(1)-0.5 xvec(end)+0.5],[0 0],'k','LineWidth',1.1)
        axis([xlim zlim])
        plot([xvec(ig)' xvec(ig)']'+0.5,[zvec(1).*ones(size(xvec(ig)')) zvec(end).*ones(size(xvec(ig)'))]','r','LineWidth',1.1)
        plot(xvec(indxn),zeros(length(indxn),1),'w*','MarkerSize',8,'LineWidth',1.1)
        n = 1;
        for i = 1:length(indx)
            temp = selmat(i,:);
            it = find(temp);
            if ~isempty(it)
                for j = 1:length(it)
                    pl(n) = plot(i,zvec(it(j)),ss(j),'Color','w','MarkerFaceColor',cs(j,:),'MarkerSize',6);
                    n = n+1;
                end
            end
        end
        mpl = find(mohovec);
        if ~isempty(mpl)
            plm = plot(mpl,mohovec(mpl),sm,'Color','w','MarkerFaceColor',cm,'MarkerSize',6);
        end
        hold off
        opt = 0;
        continue
    end

    if strcmp(anyinput,'plotallb') || strcmp(anyinput,'pab')
        if ~exist('figmain','var')
            figmain = figure('units','normalized','outerposition',[0 0.05 1 0.95]);
        else
            figure(figmain);
        end
        hold on
        imagesc(xvec,zvec,ampa2')
        colormap(jet);
        caxis([min(ampa2(:)) max(ampa2(:))])
        set(gca,'Ydir','reverse')
        plot([xvec(1)-0.5 xvec(end)+0.5],[0 0],'k','LineWidth',1.1)
        axis([xlim zlim])
        plot([xvec(ig)' xvec(ig)']'+0.5,[zvec(1).*ones(size(xvec(ig)')) zvec(end).*ones(size(xvec(ig)'))]','r','LineWidth',1.1)
        plot(xvec(indxn),zeros(length(indxn),1),'w*','MarkerSize',8,'LineWidth',1.1)
        n = 1;
        for i = 1:length(indx)
            temp = selmat(i,:);
            it = find(temp);
            if ~isempty(it)
                for j = 1:length(it)
                    pl(n) = plot(i,zvec(it(j)),ss(j),'Color','w','MarkerFaceColor',cs(j,:),'MarkerSize',6);
                    n = n+1;
                end
            end
        end
        mpl = find(mohovec);
        if ~isempty(mpl)
            plm = plot(mpl,mohovec(mpl),sm,'Color','w','MarkerFaceColor',cm,'MarkerSize',6);
        end
        hold off
        opt = 0;
        continue
    end

    if strcmp(anyinput,'plotamp1') || strcmp(anyinput,'p1')
        if ~exist('figmain','var')
            figmain = figure('units','normalized','outerposition',[0 0.05 1 0.95]);
        else
            figure(figmain);
        end
        hold on
        imagesc(xvec,zvec,amp1')
        set(gca,'Ydir','reverse')
        colormap(customColorMap);
        caxis([-max(abs(amp1(:))) max(abs(amp1(:)))])
        plot([xvec(1)-0.5 xvec(end)+0.5],[0 0],'k','LineWidth',1.1)
        axis([xlim zlim])
        plot([xvec(ig)' xvec(ig)']'+0.5,[zvec(1).*ones(size(xvec(ig)')) zvec(end).*ones(size(xvec(ig)'))]','g','LineWidth',1.1)
        plot(xvec(indxn),zeros(length(indxn),1),'g*','MarkerSize',8,'LineWidth',1.1)
        n = 1;
        for i = 1:length(indx)
            temp = selmat(i,:);
            it = find(temp);
            if ~isempty(it)
                for j = 1:length(it)
                    pl(n) = plot(i,zvec(it(j)),ss(j),'Color','k','MarkerFaceColor',cs(j,:),'MarkerSize',6);
                    n = n+1;
                end
            end
        end
        mpl = find(mohovec);
        if ~isempty(mpl)
            plm = plot(mpl,mohovec(mpl),sm,'Color','k','MarkerFaceColor',cm,'MarkerSize',6);
        end
        hold off
        opt = 1;
        continue
    end

    if strcmp(anyinput,'plotamp2') || strcmp(anyinput,'p2')
        if ~exist('figmain','var')
            figmain = figure('units','normalized','outerposition',[0 0.05 1 0.95]);
        else
            figure(figmain);
        end
        hold on
        imagesc(xvec,zvec,amp2')
        set(gca,'Ydir','reverse')
        colormap(customColorMap);
        caxis([-max(abs(amp2(:))) max(abs(amp2(:)))])
        plot([xvec(1)-0.5 xvec(end)+0.5],[0 0],'k','LineWidth',1.1)
        axis([xlim zlim])
        plot([xvec(ig)' xvec(ig)']'+0.5,[zvec(1).*ones(size(xvec(ig)')) zvec(end).*ones(size(xvec(ig)'))]','g','LineWidth',1.1)
        plot(xvec(indxn),zeros(length(indxn),1),'g*','MarkerSize',8,'LineWidth',1.1)
        n = 1;
        for i = 1:length(indx)
            temp = selmat(i,:);
            it = find(temp);
            if ~isempty(it)
                for j = 1:length(it)
                    pl(n) = plot(i,zvec(it(j)),ss(j),'Color','k','MarkerFaceColor',cs(j,:),'MarkerSize',6);
                    n = n+1;
                end
            end
        end
        mpl = find(mohovec);
        if ~isempty(mpl)
            plm = plot(mpl,mohovec(mpl),sm,'Color','k','MarkerFaceColor',cm,'MarkerSize',6);
        end
        hold off
        opt = 2;
        continue
    end

    if strcmp(anyinput,'plotamp3') || strcmp(anyinput,'p3')
        if ~exist('figmain','var')
            figmain = figure('units','normalized','outerposition',[0 0.05 1 0.95]);
        else
            figure(figmain);
        end
        hold on
        imagesc(xvec,zvec,amp3')
        set(gca,'Ydir','reverse')
        colormap(customColorMap);
        caxis([-max(abs(amp3(:))) max(abs(amp3(:)))])
        plot([xvec(1)-0.5 xvec(end)+0.5],[0 0],'k','LineWidth',1.1)
        axis([xlim zlim])
        plot([xvec(ig)' xvec(ig)']'+0.5,[zvec(1).*ones(size(xvec(ig)')) zvec(end).*ones(size(xvec(ig)'))]','g','LineWidth',1.1)
        plot(xvec(indxn),zeros(length(indxn),1),'g*','MarkerSize',8,'LineWidth',1.1)
        n = 1;
        for i = 1:length(indx)
            temp = selmat(i,:);
            it = find(temp);
            if ~isempty(it)
                for j = 1:length(it)
                    pl(n) = plot(i,zvec(it(j)),ss(j),'Color','k','MarkerFaceColor',cs(j,:),'MarkerSize',6);
                    n = n+1;
                end
            end
        end
        mpl = find(mohovec);
        if ~isempty(mpl)
            plm = plot(mpl,mohovec(mpl),sm,'Color','k','MarkerFaceColor',cm,'MarkerSize',6);
        end
        hold off
        opt = 3;
        continue
    end

    if strcmp(anyinput,'plotamp4') || strcmp(anyinput,'p4')
        if ~exist('figmain','var')
            figmain = figure('units','normalized','outerposition',[0 0.05 1 0.95]);
        else
            figure(figmain);
        end
        hold on
        imagesc(xvec,zvec,amp4')
        set(gca,'Ydir','reverse')
        colormap(customColorMap);
        caxis([-max(abs(amp4(:))) max(abs(amp4(:)))])
        plot([xvec(1)-0.5 xvec(end)+0.5],[0 0],'k','LineWidth',1.1)
        axis([xlim zlim])
        plot([xvec(ig)' xvec(ig)']'+0.5,[zvec(1).*ones(size(xvec(ig)')) zvec(end).*ones(size(xvec(ig)'))]','g','LineWidth',1.1)
        plot(xvec(indxn),zeros(length(indxn),1),'g*','MarkerSize',8,'LineWidth',1.1)
        n = 1;
        for i = 1:length(indx)
            temp = selmat(i,:);
            it = find(temp);
            if ~isempty(it)
                for j = 1:length(it)
                    pl(n) = plot(i,zvec(it(j)),ss(j),'Color','k','MarkerFaceColor',cs(j,:),'MarkerSize',6);
                    n = n+1;
                end
            end
        end
        mpl = find(mohovec);
        if ~isempty(mpl)
            plm = plot(mpl,mohovec(mpl),sm,'Color','k','MarkerFaceColor',cm,'MarkerSize',6);
        end
        hold off
        opt = 4;
        continue
    end

    if strcmp(anyinput,'plotamp5') || strcmp(anyinput,'p5')
        if ~exist('figmain','var')
            figmain = figure('units','normalized','outerposition',[0 0.05 1 0.95]);
        else
            figure(figmain);
        end
        hold on
        imagesc(xvec,zvec,amp5')
        set(gca,'Ydir','reverse')
        colormap(customColorMap);
        caxis([-max(abs(amp5(:))) max(abs(amp5(:)))])
        plot([xvec(1)-0.5 xvec(end)+0.5],[0 0],'k','LineWidth',1.1)
        axis([xlim zlim])
        plot([xvec(ig)' xvec(ig)']'+0.5,[zvec(1).*ones(size(xvec(ig)')) zvec(end).*ones(size(xvec(ig)'))]','g','LineWidth',1.1)
        plot(xvec(indxn),zeros(length(indxn),1),'g*','MarkerSize',8,'LineWidth',1.1)
        n = 1;
        for i = 1:length(indx)
            temp = selmat(i,:);
            it = find(temp);
            if ~isempty(it)
                for j = 1:length(it)
                    pl(n) = plot(i,zvec(it(j)),ss(j),'Color','k','MarkerFaceColor',cs(j,:),'MarkerSize',6);
                    n = n+1;
                end
            end
        end
        mpl = find(mohovec);
        if ~isempty(mpl)
            plm = plot(mpl,mohovec(mpl),sm,'Color','k','MarkerFaceColor',cm,'MarkerSize',6);
        end
        hold off
        opt = 5;
        continue
    end
    
    if strcmp(anyinput,'plotraw1') || strcmp(anyinput,'pr1')
        if ~exist('figmain','var')
            figmain = figure('units','normalized','outerposition',[0 0.05 1 0.95]);
        else
            figure(figmain);
        end
        hold on
        imagesc(xvec,zvec,amp1(iraw,:)')
        set(gca,'Ydir','reverse')
        colormap(customColorMap);
        caxis([-max(abs(amp1(:))) max(abs(amp1(:)))])
        plot([xvec(1)-0.5 xvec(end)+0.5],[0 0],'k','LineWidth',1.1)
        axis([xlim zlim])
%         plot([xvec(ig)' xvec(ig)']'+0.5,[zvec(1).*ones(size(xvec(ig)')) zvec(end).*ones(size(xvec(ig)'))]','g','LineWidth',1.1)
%         plot(xvec(indxn),zeros(length(indxn),1),'g*','MarkerSize',8,'LineWidth',1.1)
        n = 1;
        for i = 1:length(indx)
            temp = selmat(iraw(i),:);
            it = find(temp);
            if ~isempty(it)
                for j = 1:length(it)
                    pl(n) = plot(i,zvec(it(j)),ss(j),'Color','k','MarkerFaceColor',cs(j,:),'MarkerSize',6);
                    n = n+1;
                end
            end
        end
        mpl = find(mohovec(iraw));
        if ~isempty(mpl)
            plm = plot(mpl,mohovec(iraw(mpl)),sm,'Color','k','MarkerFaceColor',cm,'MarkerSize',6);
        end
        hold off
        opt = -2;
        continue
    end

    if strcmp(anyinput,'plotraw2') || strcmp(anyinput,'pr2')
        if ~exist('figmain','var')
            figmain = figure('units','normalized','outerposition',[0 0.05 1 0.95]);
        else
            figure(figmain);
        end
        hold on
        imagesc(xvec,zvec,amp2(iraw,:)')
        set(gca,'Ydir','reverse')
        colormap(customColorMap);
        caxis([-max(abs(amp2(:))) max(abs(amp2(:)))])
        plot([xvec(1)-0.5 xvec(end)+0.5],[0 0],'k','LineWidth',1.1)
        axis([xlim zlim])
%         plot([xvec(ig)' xvec(ig)']'+0.5,[zvec(1).*ones(size(xvec(ig)')) zvec(end).*ones(size(xvec(ig)'))]','g','LineWidth',1.1)
%         plot(xvec(indxn),zeros(length(indxn),1),'g*','MarkerSize',8,'LineWidth',1.1)
        n = 1;
        for i = 1:length(indx)
            temp = selmat(iraw(i),:);
            it = find(temp);
            if ~isempty(it)
                for j = 1:length(it)
                    pl(n) = plot(i,zvec(it(j)),ss(j),'Color','k','MarkerFaceColor',cs(j,:),'MarkerSize',6);
                    n = n+1;
                end
            end
        end
        mpl = find(mohovec(iraw));
        if ~isempty(mpl)
            plm = plot(mpl,mohovec(iraw(mpl)),sm,'Color','k','MarkerFaceColor',cm,'MarkerSize',6);
        end
        hold off
        opt = -3;
        continue
    end

    if strcmp(anyinput,'plotraw3') || strcmp(anyinput,'pr3')
        if ~exist('figmain','var')
            figmain = figure('units','normalized','outerposition',[0 0.05 1 0.95]);
        else
            figure(figmain);
        end
        hold on
        imagesc(xvec,zvec,amp3(iraw,:)')
        set(gca,'Ydir','reverse')
        colormap(customColorMap);
        caxis([-max(abs(amp3(:))) max(abs(amp3(:)))])
        plot([xvec(1)-0.5 xvec(end)+0.5],[0 0],'k','LineWidth',1.1)
        axis([xlim zlim])
%         plot([xvec(ig)' xvec(ig)']'+0.5,[zvec(1).*ones(size(xvec(ig)')) zvec(end).*ones(size(xvec(ig)'))]','g','LineWidth',1.1)
%         plot(xvec(indxn),zeros(length(indxn),1),'g*','MarkerSize',8,'LineWidth',1.1)
        n = 1;
        for i = 1:length(indx)
            temp = selmat(iraw(i),:);
            it = find(temp);
            if ~isempty(it)
                for j = 1:length(it)
                    pl(n) = plot(i,zvec(it(j)),ss(j),'Color','k','MarkerFaceColor',cs(j,:),'MarkerSize',6);
                    n = n+1;
                end
            end
        end
        mpl = find(mohovec(iraw));
        if ~isempty(mpl)
            plm = plot(mpl,mohovec(iraw(mpl)),sm,'Color','k','MarkerFaceColor',cm,'MarkerSize',6);
        end
        hold off
        opt = -4;
        continue
    end

    if strcmp(anyinput,'plotraw4') || strcmp(anyinput,'pr4')
        if ~exist('figmain','var')
            figmain = figure('units','normalized','outerposition',[0 0.05 1 0.95]);
        else
            figure(figmain);
        end
        hold on
        imagesc(xvec,zvec,amp4(iraw,:)')
        set(gca,'Ydir','reverse')
        colormap(customColorMap);
        caxis([-max(abs(amp4(:))) max(abs(amp4(:)))])
        plot([xvec(1)-0.5 xvec(end)+0.5],[0 0],'k','LineWidth',1.1)
        axis([xlim zlim])
%         plot([xvec(ig)' xvec(ig)']'+0.5,[zvec(1).*ones(size(xvec(ig)')) zvec(end).*ones(size(xvec(ig)'))]','g','LineWidth',1.1)
%         plot(xvec(indxn),zeros(length(indxn),1),'g*','MarkerSize',8,'LineWidth',1.1)
        n = 1;
        for i = 1:length(indx)
            temp = selmat(iraw(i),:);
            it = find(temp);
            if ~isempty(it)
                for j = 1:length(it)
                    pl(n) = plot(i,zvec(it(j)),ss(j),'Color','k','MarkerFaceColor',cs(j,:),'MarkerSize',6);
                    n = n+1;
                end
            end
        end
        mpl = find(mohovec(iraw));
        if ~isempty(mpl)
            plm = plot(mpl,mohovec(iraw(mpl)),sm,'Color','k','MarkerFaceColor',cm,'MarkerSize',6);
        end
        hold off
        opt = -5;
        continue
    end

    if strcmp(anyinput,'plotraw5') || strcmp(anyinput,'pr5')
        if ~exist('figmain','var')
            figmain = figure('units','normalized','outerposition',[0 0.05 1 0.95]);
        else
            figure(figmain);
        end
        hold on
        imagesc(xvec,zvec,amp5(iraw,:)')
        set(gca,'Ydir','reverse')
        colormap(customColorMap);
        caxis([-max(abs(amp5(:))) max(abs(amp5(:)))])
        plot([xvec(1)-0.5 xvec(end)+0.5],[0 0],'k','LineWidth',1.1)
        axis([xlim zlim])
%         plot([xvec(ig)' xvec(ig)']'+0.5,[zvec(1).*ones(size(xvec(ig)')) zvec(end).*ones(size(xvec(ig)'))]','g','LineWidth',1.1)
%         plot(xvec(indxn),zeros(length(indxn),1),'g*','MarkerSize',8,'LineWidth',1.1)
        n = 1;
        for i = 1:length(indx)
            temp = selmat(iraw(i),:);
            it = find(temp);
            if ~isempty(it)
                for j = 1:length(it)
                    pl(n) = plot(i,zvec(it(j)),ss(j),'Color','k','MarkerFaceColor',cs(j,:),'MarkerSize',6);
                    n = n+1;
                end
            end
        end
        mpl = find(mohovec(iraw));
        if ~isempty(mpl)
            plm = plot(mpl,mohovec(iraw(mpl)),sm,'Color','k','MarkerFaceColor',cm,'MarkerSize',6);
        end
        hold off
        opt = -6;
        continue
    end
    
    if strcmp(anyinput,'zoom') || strcmp(anyinput,'z')
        xlimt = input('Provie new limits for x (two entries): ','s');
        zlimt = input('Provide new limits for z (two entries): ','s');
        I = strfind(xlimt,' ');
        J = strfind(zlimt,' ');
        if isempty(I) || isempty(J)
            disp('more than two entries in the format: lim1 lim2')
            continue
        end
        xlim(1) = str2double(xlimt(1:I-1));
        xlim(2) = str2double(xlimt(I+1:end));
        zlim(1) = str2double(zlimt(1:J-1));
        zlim(2) = str2double(zlimt(J+1:end));
        if round(xlim(1)-0.1) == xlim(1)
            xlim(1) = xlim(1)-0.5;
        end
        if round(xlim(2)+0.1) == xlim(2)
            xlim(2) = xlim(2)+0.5;
        end
        figure(figmain);
        axis([xlim zlim])
        continue
    end

    if strcmp(anyinput,'depthzoom') || strcmp(anyinput,'dz')
        zlimt = input('Provide new limits for z (two entries): ','s');
        J = strfind(zlimt,' ');
        if isempty(J)
            disp('more than two entries in the format: lim1 lim2')
            continue
        end
        zlim(1) = str2double(zlimt(1:J-1));
        zlim(2) = str2double(zlimt(J+1:end));
        figure(figmain);
        axis([xlim zlim])
        continue
    end

    if strcmp(anyinput,'reset') || strcmp(anyinput,'r')
        zlim = [min(zvec) max(zvec)];
        xlim = [min(xvec)-0.5 max(xvec)+0.5];
        figure(figmain);
        axis([xlim zlim])
        continue
    end

    if strcmp(anyinput,'guizoom') || strcmp(anyinput,'gz')
        figure(figmain);
        hold on
        [xlim(1),zlim(1)] = ginput(1);
        pl2 = plot(xlim(1),zlim(1),'g+','LineWidth',1.1);
        [xlim(2),zlim(2)] = ginput(1);
        xlim = sort(xlim);
        zlim = sort(zlim);
        axis([xlim zlim])
        delete(pl2);
        hold off
        continue
    end

    if strcmp(anyinput,'group') || strcmp(anyinput,'g')
        disp(['There are ' num2str(max(groupsn)) ' groups.'])
        ign = input('Which group should be shown? ');
        temp = indx(groupsn==ign);
        xlim = [min(temp)-0.5 max(temp)+0.5];
        figure(figmain);
        axis([xlim zlim])
        continue
    end

    if strcmp(anyinput,'setmoho') || strcmp(anyinput,'m')
        exitflag = 0;
        nn = 1;
        if opt < 0
            mohovecnew = mohovec(iraw);
        else
            mohovecnew = mohovec;
        end
        figure(figmain);
        hold on
        while ~exitflag
            [xtemp,ztemp,button] = ginput(1);
            if button == 1
                xn(nn) = round(xtemp);
                [~,I] = min((zvec-ztemp).^2);
                zn(nn) = zvec(I);
                if opt < 0
                    temp = mohovec(iraw(xn(nn)));
                else
                    temp = mohovec(xn(nn));
                end
                if temp==0
                    pl2(nn) = plot(xn(nn),zn(nn),'g+');
                    mohovecnew(xn(nn)) = zn(nn);
                else
                    pl2(nn) = plot(xn(nn),mohovecnew(xn(nn)),'kx');
                    mohovecnew(xn(nn)) = zn(nn);
                end
                nn = nn+1;
            elseif button == 3
                exitflag = 1;
            end
        end
        acceptflag = input('Accept new layer assignments? 1=yes,0=no');
        if ~acceptflag
            delete(pl2)
        else
            delete(pl2)
            if exist('plm','var')
            delete(plm)
            end
            mohovec = mohovecnew;
            mpl = find(mohovec);
            if opt == 0 || opt == 1
                if ~isempty(mpl)
                    plm = plot(mpl,mohovec(mpl),sm,'Color','w','MarkerFaceColor',cm,'MarkerSize',6);
                end
            else
                if ~isempty(mpl)
                    plm = plot(mpl,mohovec(mpl),sm,'Color','k','MarkerFaceColor',cm,'MarkerSize',6);
                end
            end
            if opt < 0
                mohovec = mohovec(jjn);
            end
        end
        hold off
        continue
    end

    if strcmp(anyinput,'setlayers') || strcmp(anyinput,'s')
        exitflag = 0;
        nn = 1;
        if opt < 0
            selmatnew = selmat(iraw,:);
        else
            selmatnew = selmat;
        end
        figure(figmain);
        hold on
        while ~exitflag
            [xtemp,ztemp,button] = ginput(1);
            if button == 1
                xn(nn) = round(xtemp);
                [~,I] = min((zvec-ztemp).^2);
                zn(nn) = zvec(I);
                iz = find(zvec==zn(nn),1);
                if opt < 0
                    temp = selmat(iraw(xn(nn)),:);
                else
                    temp = selmat(xn(nn),:);
                end
                izt = find(temp(iz-5:min(length(zvec),iz+5)));
                if isempty(izt)
                    pl2(nn) = plot(xn(nn),zn(nn),'g+');
                    selmatnew(xn(nn),iz) = 1;
                else
                    izt = izt(1);
                    pl2(nn) = plot(xn(nn),zvec(iz-5+izt-1),'kx');
                    selmatnew(xn(nn),iz-5+izt-1) = 0;
                end
                nn = nn+1;
            elseif button == 3
                exitflag = 1;
            end
        end
        acceptflag = input('Accept new layer assignments? 1=yes,0=no');
        if ~acceptflag
            delete(pl2)
        else
            delete(pl2)
            if exist('pl','var')
            delete(pl)
            end
            selmat = selmatnew;
            if opt == 0 || opt == -1
                n = 1;
                for i = 1:length(indx)
                    temp = selmat(i,:);
                    it = find(temp);
                    if ~isempty(it)
                        for j = 1:length(it)
                            pl(n) = plot(i,zvec(it(j)),ss(j),'Color','w','MarkerFaceColor',cs(j,:),'MarkerSize',6);
                            n = n+1;
                        end
                    end
                end
            else
                n = 1;
                for i = 1:length(indx)
                    temp = selmat(i,:);
                    it = find(temp);
                    if ~isempty(it)
                        for j = 1:length(it)
                            pl(n) = plot(i,zvec(it(j)),ss(j),'Color','k','MarkerFaceColor',cs(j,:),'MarkerSize',6);
                            n = n+1;
                        end
                    end
                end
            end
            if opt < 0
                selmat = selmat(jjn,:);
            end
        end
        hold off
        continue
    end
    
    if strcmp(anyinput,'help')
        disp('List of commands:');
        disp('    - quit          : close program')
        disp('    - plotraw/pr       : plot figure with squared and summed harmonics (original order)')
        disp('    - plotall/pa       : plot figure with squared and summed harmonics (cluster groups)')
        disp('    - plotamp1/p1      : plot figure with constant term')
        disp('    - plotamp2/p2      : plot figure with cos(baz) term')
        disp('    - plotamp3/p3      : plot figure with sin(baz) term')
        disp('    - plotamp4/p4      : plot figure with cos(2*baz) term')
        disp('    - plotamp5/p5      : plot figure with sin(2*baz) term')
        disp('    - zoom/z           : zoom into the figure by entering x and z limits')
        disp('    - guizoom/gz       : zoom into the figure by graphically setting boundary box')
        disp('    - group/g          : show only specific group from cluster analysis')
        disp('    - setlayer/s       : assign layers or delete assignments')
        disp('    - reset/r          : reset zoom')
        continue
    end

    try 
        eval(anyinput)
    catch
        disp('No matlab or PI command')
    end
end

nn = 1;
nlay = [];
for i = 1:length(indx)
    temp = selmat(indx(i),:);
    if mohovec(indx(i)) == 0
        continue
    end
    if sum(temp) == 0 
        continue
    else
        stats{nn} = stat{indx(i)};
        mohoz(nn) = mohovec(indx(i));
        ix = find(temp);
        layt{nn} = zvec(ix);
        if isempty(nlay)
            nlay = length(ix);
        else
            nlay = max(nlay,length(ix));
        end
        nn = nn+1;
    end
end
lays = zeros(length(stats),nlay);
for i = 1:length(stats)
    lays(i,1:length(layt{i})) = layt{i};
end
filename = 'LayerModel.txt';
fid = fopen([filedir '\' filename],'wt');
% format = ['%s,' repmat('%f,',1,nlay) '%f,%f,%f\n'];
% for i = 1:length(stats)
%     fprintf(fid,format,stats{i},mohoz(i),lays(i,:),groupsn(i),~isempty(find(indxn==indx(i),1)));
% end
format = ['%s,' repmat('%f,',1,nlay) '%f\n'];
for i = 1:length(stats)
    fprintf(fid,format,stats{i},mohoz(i),lays(i,:));
end
fclose(fid);
filenamefig = 'LayerModel';
print(figmain,'-r600','-djpeg',[filedir '/' filenamefig '.jpg'])
print(figmain,'-r600','-dpdf','-bestfit',[filedir '/' filenamefig '.pdf'])