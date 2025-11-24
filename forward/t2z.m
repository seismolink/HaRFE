function z_out = t2z(ttarget,ldepth,velp,vels,slow)
% Map depth to time   
    nlayers = length(ldepth)-1;
    dt = zeros(size(ldepth));
    zout = zeros(size(ttarget));


    ddepth(1) = ldepth(1);
    coss = sqrt(1-vels(1)*vels(1)*slow*slow);
    cosp = sqrt(1-velp(1)*velp(1)*slow*slow);
    dt(1) = ddepth(1)*(coss/vels(1)-cosp/velp(1));
    lt(1) = dt(1);
    for ii = 1:nlayers
        ddepth(ii+1) = ldepth(ii+1)-ldepth(ii);
        coss = sqrt(1-vels(ii+1)*vels(ii+1)*slow*slow);
        cosp = sqrt(1-velp(ii+1)*velp(ii+1)*slow*slow);
        dt(ii+1) = ddepth(ii+1)*(coss/vels(ii+1)-cosp/velp(ii+1));
        lt(ii+1) = sum(dt);
    end

    for i = 1:length(ttarget)
    zshft = 0; %tshft = 0;
    mlay = 1;
    rest = ttarget(i);
    for ii = 1:length(dt)
        t = lt(ii);
        if ttarget(i) > t %ztarget > depth
            mlay = mlay+1;
            rest = ttarget(i) - t;
        end
    end
    if mlay > 1
        for k = 1:mlay-1
            coss = sqrt(1-vels(k)*vels(k)*slow*slow);
            cosp = sqrt(1-velp(k)*velp(k)*slow*slow);
            zshft = zshft+ dt(k)/(coss/vels(k)-cosp/velp(k));
        end
    end

    coss = sqrt(1-vels(mlay)*vels(mlay)*slow*slow);
    cosp = sqrt(1-velp(mlay)*velp(mlay)*slow*slow);
    zshft = zshft+ rest/(coss/vels(mlay)-cosp/velp(mlay));
    z_out(i) = zshft;
    end
end