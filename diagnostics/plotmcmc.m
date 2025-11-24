clear all
close all
A = dir('D:\MARSn\MARSa\MCMCresult*.mat');
id1 = 18;
id2 = 21;
% A = dir('C:\Users\frede\OneDrive\Dokumente\WorkDocs\Documents\YalePostDoc\HarmDecMCMC\Inversions\Bourke\BCIP_v1\MCMCresult*.mat');
% A = dir('C:\Users\frede\OneDrive\Dokumente\WorkDocs\Documents\YalePostDoc\HarmDecMCMC\Inversions\SeisConn8\CS07\MCMCresult*.mat');
% id1 = 22;
% id2 = 25;
figure; hold on; 

for i = 1:length(A)
    load([A(i).folder,'/',A(i).name]);
    plot(MC.post.full.stat(:,id1),MC.post.full.stat(:,id2),'k-','LineWidth',0.01); 
end
for i = 1:length(A)
    load([A(i).folder,'/',A(i).name]);
    scatter(MC.post.full.stat(:,id1),MC.post.full.stat(:,id2),1.5*(log(1./(1:length(MC.post.full.stat(:,id1))))+10),exp(MC.post.full.L1),'filled')
end
colorbar