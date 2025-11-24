function [tt,amp,p,evecin_list,amp_list,tt_list,bailout]=raysum(thick,rho,isoflag,aa,ar_list,rot,baz,slow,phase,nseg,nlay,amp_in,p,evecin_list,amp_list,tt_list,fchg)
% Main arrival-calculation routine. See raysum2.m
% Update -- keeps complex matrices as long as possible, for
% accuracy. p and amplitude are kept real -- complex values
% make no sense in ray theory.

% Interface variables:
 
bailout = false;

if(fchg == 0)
    % Initialize variables, get incident slowness vector
    % Assume half-space is isotropic -- save some grief.
    tt_list(1)=0;
    amp_list(1)=amp_in;
    
    p(1,1)=-slow.*cos(baz);
    p(2,1)=-slow.*sin(baz);
    
    [evalbot,evecbot] = Hisotroc(aa(:,:,:,:,nlay),rho(nlay),p(1,1),p(2,1));
    phase1=fix(rem(phase(1,2)+2,6)+1);
    p(3,1)=real(evalbot(phase1));

    % Store active eigenvector for later sign check.
    evecin_list(:,1) = evecbot(1:3,phase1);
end

for seg=max((fchg-1),1):(nseg-1)

% 1 is incident, 2 is transmitted/refelcted
lay1=phase(seg,1);
lay2=phase(seg+1,1);
% First and second phase, converted to eigenvector index
phase1=fix(rem(phase(seg,2)+2,6)+1);
phase2=fix(rem(phase(seg+1,2)+2,6)+1);
% Distinguish upper and lower layers.
laytop=min(lay1,lay2);
laybot=max(lay1,lay2);
% Upflag: true_ml if the incident wave is upgoing.
upflag=phase1 > 3;
% Rflag: true_ml for a reflection.
rflag=lay1 == lay2;

if(rflag)
    if(upflag)
        laytop=laybot-1;
    else
        laybot=laytop+1;
    end
end

% Fsflag: free-surface reflection
fsflag=laytop == 0;
if(fsflag &&(~ rflag))
    % disp('can''t have free surface transmission')
end

% laytop and laybot now contain correct layers for R/T calculation.
% Find correct rotator, determining active interface:
if(upflag)
    rnum=lay1;
else
    rnum=lay1+1;
end
rt = rot(:,:,rnum)';

% Project slowness onto interface
p1r = rt*p(:,seg);

% If the slowness vector is incorrectly oriented wrt the interface,
% we have a trapped phase.
if((upflag &&(p1r(3) > 0)) || ((~ upflag) &&(p1r(3) < 0)))
    % disp('Ray does not intersect interface. No. 1')
    amp(1)=0;
    amp(2)=0;
    amp(3)=0;
    tt=0;
    tt_list(seg+1)=0;
    amp_list(seg+1)=0;
    bailout = true;
    return
end

% Find evals/evecs. This is where it gets hairy...

%    Lower eigenvalues/vectors
if(isoflag(laybot))
    [evalbot,evecbot] = Hisotroc(aa(:,:,:,:,laybot),rho(laybot),p1r(1),p1r(2));
else
    %           Upper interface of laybot
    [evalbot,evecbot] = Hanisotroc(ar_list(:,:,:,:,laybot,2),rho(laybot),p1r(1),p1r(2));
end

%    Handle degeneracy (lower)
if (abs((evalbot(2)-evalbot(3))./evalbot(2)) <= eps)
    [evecbot] = Hrot_evec(evecbot,rot(:,:,rnum));
end

%     Upper eigenvalues/vectors. Don't bother if we're at the
%     free surface
if (~ fsflag)
if (isoflag(laytop))
    [evaltop,evectop] = Hisotroc(aa(:,:,:,:,laytop),rho(laytop),p1r(1),p1r(2));
else
%           Lower interface of laytop
    [evaltop,evectop] = Hanisotroc(ar_list(:,:,:,:,laytop,1),rho(laytop),p1r(1),p1r(2));
end
%              Handle degeneracy (upper)
if ((evaltop(2) ~= 0.) && (abs((evaltop(2)-evaltop(3))./evaltop(2)) <= eps))
    [evectop] = Hrot_evec(evectop,rot(:,:,rnum));
end
end

% Distribute evals/evecs between first and second propagation legs:

if(rflag)
%          Reflection
    if(upflag)
    %            Upgoing wave reflects off top layer -- incl. free-surf.
    %            evec1=evecbot; eval1=evalbot; evec2=evec1; eval2=eval1;
        eval1 = evalbot;
        evec1 = evecbot;
        eval2 = eval1;
        evec2 = evec1;
    else
    %            Downgoing wave reflects off bottom layer.
    %            evec1=evectop; eval1=evaltop; evec2=evec1; eval2=eval1;
        eval1 = evaltop;
        evec1 = evectop;
        eval2 = eval1;
        evec2 = evec1;
    end
else
    %          Transmission
    if(upflag)
    %            Upgoing wave enters top layer
    %            evec1=evecbot; eval1=evalbot; evec2=evectop; eval2=evaltop;
        eval1 = evalbot;
        evec1 = evecbot;
        eval2 = evaltop;
        evec2 = evectop;
    else
    %            Downgoing wave enters bottom layer
    %            evec1=evectop; eval1=evaltop; evec2=evecbot; eval2=evalbot;
        eval1 = evaltop;
        evec1 = evectop;
        eval2 = evalbot;
        evec2 = evecbot;
    end
end

%  If the eigenvalue is zero, we didn't get a transmission -- bail out
if(abs(real(eval2(phase2))) < eps^2)
    % disp('Phase not transmitted')
    amp(1)=0;
    amp(2)=0;
    amp(3)=0;
    tt=0;
    tt_list(seg+1)=0;
    amp_list(seg+1)=0;
    bailout = true;
    return
end


% Consistency check: are the incident phase and sign correct?
% Consistency may regained by altering phase1 or using a sign
% multiplier (mult)
evecin_r = rt*evecin_list(:,seg);
[phase1,mult_ml,errflag] = Hevec_check(evec1,evecin_r,phase1);
if(errflag)
%     disp('Error')
end

% Find new slowness vector. By Snell's law, interface-parallel
% components carry over. Eigenvalue provides remaining component.
p2r(1)=p1r(1);
p2r(2)=p1r(2);
p2r(3)=real(eval2(phase2));
p(:,seg+1) = rot(:,:,rnum)*p2r';

% Now we're ready to calculate the amplitude. Two cases to consider.
if (fsflag)
%          Free-surface multiple
%            Partition N: Nu=evecbot(4:6,4:6), Nd=evecbot(4:6,1:3)
    nu = evecbot(4:6,4:6);
    nd = evecbot(4:6,1:3);
% Free-surface reflection matrix: MM=-Nd^(-1)*Nu
    mm = -nd\nu;
else
%          Layer interaction
%            scattering matrix QQ=evecbot^(-1)*evectop
    qq = evecbot\evectop;
%            Pull out reflection or transmission matrix MM
    if (rflag)
        if (upflag)
        %                Ru = Q(1:3,4:6)*inv(Q(4:6,4:6))
            mm = qq(1:3,4:6)/qq(4:6,4:6);
        else
        %                Rd = -inv(Q(4:6,4:6))*Q(4:6,1:3)
            mm = -qq(4:6,4:6)\qq(4:6,1:3);
        end
    else
        if (upflag)
        %                Tu=inv(Q(4:6,4:6))
            mm = inv(qq(4:6,4:6));
        else
        %                Td=Q(1:3,1:3)-Q(1:3,4:6)*inv(Q(4:6,4:6))*Q(4:6,1:3)
            mm = qq(1:3,1:3)-(qq(1:3,4:6)/qq(4:6,4:6))*qq(4:6,1:3);
        end
    end
end

%       Carry amplitude across interface
amp_list(seg+1)=amp_list(seg).*real(mult_ml).*real(mm(rem(phase2-1,3)+1,rem(phase1-1,3)+1));

%       Add post-interface segment to travel time
if (thick(lay2) < 0.)
%     disp('ErroR')
end
tt_list(seg+1)=abs(p(3,seg+1)).*thick(lay2);

%       Pull out relevant eigenvector for next-segment consistency
%       check. evecin=evec_out(:,seg+1)=R*evec2(1:3,phase2)
evecin_list(:,seg+1) = rot(:,:,rnum)*evec2(1:3,phase2);
end

%      Propagation complete. Now just need to convert amplitude to
%      displacement, using the free-surface transfer matrix.
%      MAY BE BUGGY FOR ANISOTROPIC TOPMOST LAYER.

%         Check if last segment really is upgoing. Otherwise, bail out.
if (p(3,nseg) > 0.)
    % disp('Ray does not reach surface. No. 2')
    amp(1)=0;
    amp(2)=0;
    amp(3)=0;
    tt=0;
    return
end

%    Eigenvectors, again:
if (isoflag(1))
    [~,evectop] = Hisotroc(aa(:,:,:,:,1),rho(1),p(1,nseg), p(2,nseg));
else
    [~,evectop] = Hanisotroc(aa(:,:,:,:,1),rho(1),p(1,nseg), p(2,nseg));
end
phase1=fix(rem(phase(nseg,2)+2,6)+1);
[phase1,mult_ml,errflag] = Hevec_check(evectop,evecin_list(:,nseg),phase1);
if(errflag)
%     disp('error')
end

%  Get amplitude, from free-surface transfer, as follows:
%  Md=evec1(1:3,1:3),Mu=evec1(1:3,4:6),Nd=evec1(4:6,1:3),Nu=evec1(4:6,4:6)
%      cu(mod(phase1-1,3)+1)=amp_out(nseg)*mult;
%      amp=-(Mu-Md*inv(Nd)*Nu)*cu; p1r stands for cu (wavevector)
p1r(1)=0;
p1r(2)=0;
p1r(3)=0;
p1r(rem(phase1-1,3)+1)=amp_list(nseg).*mult_ml;
amp = -(evectop(1:3,4:6)-(evectop(1:3,1:3)/evectop(4:6,1:3))*evectop(4:6,4:6))*p1r;

% Assemble travel-time:
tt = sum(tt_list(1:nseg));

end