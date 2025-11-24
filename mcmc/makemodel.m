function model = makemodel(idx,b)
% Produce a model for forwardHDref from a given normalized model vector idx
% (9 values between 0 and 1) and a boundary structure with min and max
% bounds for each variable
n = size(idx);
model(1,10) = 0;
model(1,11) = 0;
for i = 1:n(1)
    if idx(i,1) == -1
        model(i,1) = model(i-1,1);
    else
        model(i,1) = idx(i,1)*(b.hmax(i)-b.hmin(i))+b.hmin(i);
    end
    if idx(i,4) == -1
        model(i,2) = model(i-1,2);
    else
        model(i,2) = idx(i,4)*(b.rhomax(i)-b.rhomin(i))+b.rhomin(i);
    end
    if idx(i,2) == -1
        model(i,3) = model(i-1,3);
    else
        model(i,3) = idx(i,2)*(b.vpmax(i)-b.vpmin(i))+b.vpmin(i);
    end
    % if idx(i,3) == -1
    %     model(i,4) = model(i-1,4);
    % else
    %     model(i,4) = model(i,3)./(idx(i,3)*(b.vsmax(i)-b.vsmin(i))+b.vsmin(i));
    % end
    if idx(i,3) == -1
        model(i,4) = model(i-1,4);
    else
        model(i,4) = idx(i,3)*(b.vsmax(i)-b.vsmin(i))+b.vsmin(i);
    end
    if idx(i,5) == -1
        if mod(randi(2,1),2)==1
            model(i,6) = model(i-1,6);
            model(i-1,6) = 0;
            flag = 1;
        else
            model(i,6) = 0;
            flag = 0;
        end
    else
        model(i,6) = idx(i,5)*(b.amax(i)-b.amin(i))+b.amin(i);
    end
    model(i,7) = model(i,6);
    model(i,5) = model(i,6)==0;
    if idx(i,6) == -1
        if flag == 1
            model(i,8) = model(i-1,8);
            model(i-1,8) = 0;
        else
            model(i,8) = 0;
        end
    else
        model(i,8) = idx(i,6)*(b.phimax(i)-b.phimin(i))+b.phimin(i);
    end
    if idx(i,7) == -1
        if flag == 1
            model(i,9) = model(i-1,9);
            model(i-1,9) = 0;
        else
            model(i,9) = 0;
        end
    else
        model(i,9) = idx(i,7)*(b.plmax(i)-b.plmin(i))+b.plmin(i);
    end
    if idx(i,8) == -1
        if flag == 1
            model(i+1,10) = model(i,10);
            model(i,10) = 0;
        else
            model(i+1,10) = 0;
        end
    else
        model(i+1,10) = idx(i,8)*(b.srmax(i)-b.srmin(i))+b.srmin(i);
    end
    if idx(i,9) == -1
        if flag == 1
            model(i+1,11) = model(i,11);
            model(i,11) = 0;
        else
            model(i+1,11) = 0;
        end
    else
        model(i+1,11) = idx(i,9)*(b.dpmax(i)-b.dpmin(i))+b.dpmin(i);
    end
    
%     model(:,1) = [idx(1)*(max(b.h)-min(b.h))+min(b.h) 0]'.*1000; %h
%     model(:,2) = [idx(2)*(max(b.rho)-min(b.rho))+min(b.rho) 3600]'; %rho
%     model(:,3) = [idx(3)*(max(b.vp)-min(b.vp))+min(b.vp) idx(3)*(max(b.vp)-min(b.vp))+min(b.vp)+idx(10)*(max(b.dvp)-min(b.dvp))+min(b.dvp)]'.*1000; %vp
%     model(:,4) = [idx(4)*(max(b.vs)-min(b.vs))+min(b.vs) (idx(3)*(max(b.vp)-min(b.vp))+min(b.vp)+idx(10)*(max(b.dvp)-min(b.dvp))+min(b.dvp))/1.8]'.*1000; %vs
%     model(:,5) = [(idx(5)*(max(b.a)-min(b.a))+min(b.a))==0 1]'; %isoflag
%     model(:,6) = [idx(5)*(max(b.a)-min(b.a))+min(b.a) 0]'; %avp
%     model(:,7) = [idx(5)*(max(b.a)-min(b.a))+min(b.a) 0]'; %avs
%     model(:,8) = [idx(6)*(max(b.phi)-min(b.phi))+min(b.phi) 0]'; %phi
%     model(:,9) = [idx(7)*(max(b.pl)-min(b.pl))+min(b.pl) 0]'; %plunge
end
model(2:end-1,1) = diff(model(1:end-1,1));
model(n(1)+1,1) = 0;
model(n(1)+1,2) = b.rhm;
if ~isfield(b,'pm')
model(n(1)+1,3) = model(n(1),3)+idx(1,10)*(b.dvp(2)-b.dvp(1))+b.dvp(1);
model(n(1)+1,4) = model(n(1)+1,3)/b.km;
else    
if b.pm == 0
model(n(1)+1,3) = model(n(1),3)+idx(1,10)*(b.dvp(2)-b.dvp(1))+b.dvp(1);
model(n(1)+1,4) = model(n(1)+1,3)/b.km;
else
% model(n(1)+1,3) = b.pm;
% model(n(1)+1,4) = model(n(1),4)+idx(1,10)*(b.dvp(2)-b.dvp(1))+b.dvp(1);
model(n(1)+1,3) = model(n(1),3)+idx(1,10)*(b.dvp(2)-b.dvp(1))+b.dvp(1);
model(n(1)+1,4) = model(n(1)+1,3)/b.km;
end
end
model(n(1)+1,5) = 1;
model(n(1)+1,6) = 0;
model(n(1)+1,7) = 0;
model(n(1)+1,8) = 0;
model(n(1)+1,9) = 0;


model(abs(model(:,6))<1,7) = 0;
model(abs(model(:,6))<1,8) = 0;
model(abs(model(:,6))<1,9) = 0;
model(abs(model(:,6))<1,5) = 1;
model(abs(model(:,6))<1,6) = 0;


model(:,1:4) = model(:,1:4)*1000;
end