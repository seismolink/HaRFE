function [travel_time,amplitude]=get_arrivals(thick,rho,isoflag,strike,dip,aa,ar_list,rot,baz,slow,sta_dx,sta_dy,phaselist,ntr,nseg,numph,nlay,amp_in)
% Get a list of arrivals, given model, geometry, and desired phases.

% Calculate thickness changes caused by shifts from baseline
dzdx = tan(dip).*sin(-strike);
dzdy = tan(dip).*sin(pi./2-strike);

dthdx = diff(dzdx);
dthdy = diff(dzdy);
dthdx(end+1) = 0;
dthdy(end+1) = 0;

travel_time = zeros(numph,ntr);
amplitude = zeros(3,numph,ntr);
amp1 = zeros(numph,ntr);
amp2 = amp1;
amp3 = amp1;
amptemp = zeros(3,1);
for itr=1:ntr
    thick_shift=thick+sta_dx(itr).*dthdx+sta_dy(itr).*dthdy;
    
    bailout = true;
    for iph=1:numph
        if(bailout)
            fchg=0;
            p = zeros(3,1);
            evecin_list = zeros(3,1);
            amp_list = zeros(3,1);
            tt_list = 0;
        else
            fchg=1;
            if(nseg(iph) == nseg(iph-1))
                while((phaselist(fchg,2,iph) == phaselist(fchg,2,iph-1))&&(phaselist(fchg,1,iph) == phaselist(fchg,1,iph-1)))
                    fchg=fix(fchg+1);
                end
            end
        end
%         [travel_time(iph,itr),amplitude(1:end,iph,itr),p,evecin_list,amp_list,tt_list,bailout]=Hraysum(thick_shift,rho,isoflag,aa,ar_list,rot,baz(itr),slow(itr),phaselist(1:end,1:end,iph),nseg(iph), nlay,amp_in,p,evecin_list,amp_list,tt_list,fchg);
        [travel_time(iph,itr),amptemp,p,evecin_list,amp_list,tt_list,bailout]=Hraysum(thick_shift,rho,isoflag,aa,ar_list,rot,baz(itr),slow(itr),phaselist(1:end,1:end,iph),nseg(iph), nlay,amp_in,p,evecin_list,amp_list,tt_list,fchg);
        amp1(iph,itr) = amptemp(1);
        amp2(iph,itr) = amptemp(2);
        amp3(iph,itr) = amptemp(3);
    end
end
amplitude(1,:,:) = amp1;
amplitude(2,:,:) = amp2;
amplitude(3,:,:) = amp3;
end