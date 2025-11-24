function t_out = z2t(ztarget,ldepth,velp,vels,slow)
% Map depth to time    
    nlayers = length(ldepth)-1;
    ddepth = zeros(size(ldepth));

    ddepth(1) = ldepth(1);
    for ii = 1:nlayers
        ddepth(ii+1) = ldepth(ii+1)-ldepth(ii);
    end

    tshft = 0;
    mlay = 1;
    resz = ztarget;
    for ii = 1:length(ldepth)
        depth = ldepth(ii);
        if ztarget > depth
            mlay = mlay+1;
            resz = ztarget - depth;
        end
    end
    if mlay > 1
        for k = 1:mlay-1
            coss = sqrt(1-vels(k)*vels(k)*slow*slow);
            cosp = sqrt(1-velp(k)*velp(k)*slow*slow);
            tshft = tshft+ ddepth(k)*(coss/vels(k)-cosp/velp(k));
        end
    end

    coss = sqrt(1-vels(mlay)*vels(mlay)*slow*slow);
    cosp = sqrt(1-velp(mlay)*velp(mlay)*slow*slow);
    tshft = tshft+ resz*(coss/vels(mlay)-cosp/velp(mlay));
    t_out = tshft;
end