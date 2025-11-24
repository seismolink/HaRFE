clear all
close all

filedir = input('Define path to results: ','s');
programdir = 'C:\Users\frede\OneDrive\Dokumente\WorkDocs\Documents\YalePostDoc\HarmDecMCMC\HDinv4\';

filename = 'LayerModel.txt';
fid = fopen([filedir '/' filename],'rt');
line = fgetl(fid);
n = 1;
while ischar(line)
    ix = strfind(line,',');
    format = ['%s ' repmat('%f ',1,length(ix)-1), '%f'];
    C = textscan(line,format,'Delimiter',',');

    statoi{n} = char(C{1,1});
    moho{n} = C{1,2};
    for i = 2:length(ix)
        if C{1,i+1} == 0
            break;
        end
        lays{n}(i-1) = C{1,i+1};
    end
    formatnew{n} = [format((1:(length(lays{n})*3))+3) ' %f'];
    % disp(formatnew{n})
    n = n+1;
    line = fgetl(fid);
end
N = length(statoi);
fclose(fid);

% fid = fopen([programdir '/mcmc/PEMC.csv'],'rt');
% C = textscan(fid,'%f %f %f %f %f','Delimiter',',');
% fclose(fid);
% depth = C{1,2};
% rh = C{1,3};
% vp = C{1,4};
% vs = C{1,5};
% vp(diff(vp)==0)

% cite: 
fid = fopen([programdir '/mcmc/App1h5.csv'],'rt');
C = textscan(fid,'%f %f %f %f %f %f %f','Delimiter',',');
fclose(fid);
depth = C{1,2};
rh = C{1,5};
vp = C{1,3};
vs = C{1,4};


% load([programdir '/mcmc/ak135.mat'],'speed');
% depth0=speed(:,1); %%convert everything to SI units
% vp0=speed(:,3);
% vs0=speed(:,2);
% rh0=speed(:,4);
% h00 = [0 depth0'];
% depth00(2) = 2;
% vp00 = [3.5 vp0'];
% vs00 = [3.5/1.73 vs0'];
% rh00 = [2.5 rh0'];
% n = 1;
% depth(n) = depth00(1);
% vp(n) = vp00(1)-0.1;
% vs(n) = vs00(1)-0.1;
% rh(n) = rh00(1)-0.1;
% n = n+1;
% for i = 1:length(depth00)-1
%     depth(n) = depth00(i+1);
%     vp(n) = vp00(i);
%     vs(n) = vs00(i);
%     rh(n) = rh00(i);
%     n = n+1;
%     depth(n) = depth00(i+1)+0.1;
%     vp(n) = vp00(i+1);
%     vs(n) = vs00(i+1);
%     rh(n) = rh00(i+1);
%     n = n+1;
% end
% dept

stats = {'CS05','CS06','CS07','CS08','CS09'};
vps = 3.5;
vss = 2.0;
rhs = 2.5;

filename2 = [filedir '/RhoVpVsInit.txt'];

rhc = interp1(depth(1:8),rh(1:8),linspace(depth(1),depth(8),100));
vpc = interp1(depth(1:8),vp(1:8),linspace(depth(1),depth(8),100));
vsc = interp1(depth(1:8),vs(1:8),linspace(depth(1),depth(8),100));
dc = interp1(depth(1:8),depth(1:8),linspace(depth(1),depth(8),100));
rhm = interp1(depth(9:end),rh(9:end),linspace(depth(1),depth(end),10000));
vpm = interp1(depth(9:end),vp(9:end),linspace(depth(1),depth(end),10000));
vsm= interp1(depth(9:end),vs(9:end),linspace(depth(1),depth(end),10000));
dm = interp1(depth(9:end),depth(9:end),linspace(depth(1),depth(end),10000));
for i = 1:length(moho)
    clear vpe vse rhe
    mt = moho{i};
    laymod = lays{i};
depoi = mt;
depest = sort([laymod mt]);
% z0 = 0:0.1:depoi+0.1;

z0 = linspace(0,depoi+0.1,100);
for j = 1:length(depest(depest<=mt))
    if j == 1
        vpe(j) = mean(vpc(z0<depest(j)));
        vse(j) = mean(vsc(z0<depest(j)));
        rhe(j) = mean(rhc(z0<depest(j)));
    else
        vpe(j) = mean(vpc(z0<depest(j)&z0>=depest(j-1)));
        vse(j) = mean(vsc(z0<depest(j)&z0>=depest(j-1)));
        rhe(j) = mean(rhc(z0<depest(j)&z0>=depest(j-1)));
    end
end
J = length(vpe);
for j = 1:length(depest(depest>mt))
    vpe(J+j) = mean(vpm(dm<depest(J+j)&dm>=depest(J+j-1)));
    vse(J+j) = mean(vsm(dm<depest(J+j)&dm>=depest(J+j-1)));
    rhe(J+j) = mean(rhm(dm<depest(J+j)&dm>=depest(J+j-1)));
end
if depest(end) < dm(1)
    vpe(end+1) = vpm(1);
    vse(end+1) = vsm(1);
    rhe(end+1) = rhm(1);
else
temp = vpm(dm>=depest(end));
vpe(end+1) = temp(2);
temp = vsm(dm>=depest(end));
vse(end+1) = temp(2);
temp = rhm(dm>=depest(end));
rhe(end+1) = temp(2);
% return
% 
% % prepare vp model
% 
% % load([programdir '/mcmc/ak135.mat'],'speed');
% % depth=speed(:,1); %%convert everything to SI units
% % vp=speed(:,3);
% % vs=speed(:,2);
% % rh=speed(:,4);
% % vp(1:2) = [-0.2 -0.1]+vp(3);
% % vs(1:2) = [-0.2 -0.1]+vs(3);
% % rh(1:2) = [-0.2 -0.1]+rh(3);
% z2 = 0:0.1:150;
% for i = 2:length(z2)
%     vpi2(i) = interp1(depth,vp,z2(i));
%     vsi2(i) = interp1(depth,vs,z2(i));
%     rhi2(i) = interp1(depth,rh,z2(i));
% end
% vpe(depest>max(z0)) = interp1(z2,vpi2,depest(depest>max(z0)));
% vse(depest>max(z0)) = interp1(z2,vsi2,depest(depest>max(z0)));
% rhe(depest>max(z0)) = interp1(z2,rhi2,depest(depest>max(z0)));
if contains(stats,statoi{i})
    vpe(1) = vps;
    vse(1) = vss;
    rhe(1) = rhs;
end
vpe = round(vpe*100/5)*5/100;
vse = round(vse*100/5)*5/100;
rhe = round(rhe*100/5)*5/100;
format = ['%s,' repmat('%f,',1,length(vpe)+length(vse)+length(rhe)), '%f\n'];
fid = fopen(filename2,'at');
fprintf(fid,format,statoi{i},rhe,vpe,vse,1);
fclose(fid);
end
end

return




preflag2 = 0;
if exist([filedir '/' filename2],'file')
    fid = fopen([filedir '/' filename2],'rt');
    line = fgetl(fid);
    n = 1;
    while ischar(line)
        ix = strfind(line,',');
        format = ['%s ' repmat('%f ',1,length(ix)-1), ' %f'];
        C = textscan(line,format,'Delimiter',',');
   	 
        statoi2{n} = char(C{1,1});
        for i = 2:length(ix)
            pars{n}(i-1) = C{1,i};
        end
        setvpflag{n} = C{1,end};
        % disp(num2str(setvpflag{n}))
	n = n+1;
        line = fgetl(fid);
    end
    preflag2 = 1;
end