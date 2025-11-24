function tout = z2t(ztarget,ldepth,velp,vels,slow)
% Map depth to time    
nlayers = length(ldepth)-1;
ddepth = zeros(size(ldepth));
tout = zeros(size(ztarget));

ddepth(1) = ldepth(1);
for ii = 1:nlayers
    ddepth(ii+1) = ldepth(ii+1)-ldepth(ii);
end

for i = 1:length(ztarget)
tshft = 0;
mlay = 0;
resz = ztarget(i);
for ii = 1:length(ldepth)
    depth = ldepth(ii);
    if ztarget(i)>depth;
        mlay = mlay+1;
        resz = ztarget(i)-depth;
    end
end

if mlay > 0
    for k = 1:mlay
        coss = sqrt(1-vels(k)*vels(k)*slow*slow);
        cosp = sqrt(1-velp(k)*velp(k)*slow*slow);
        tshft = tshft+ddepth(k)*(coss/vels(k)-cosp/velp(k));
    end
end

coss = sqrt(1-vels(mlay+1)*vels(mlay+1)*slow*slow);
cosp = sqrt(1-velp(mlay+1)*velp(mlay+1)*slow*slow);
tshft = tshft+resz*(coss/vels(mlay+1)-cosp/velp(mlay+1));
tout(i) = tshft;
end
end