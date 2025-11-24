% plot all profile stuff
clear all
close all

projfolder = 'D:\flink\HimalayaAnalysis\Inversion3\';
freq = '1.0';
addname = 'FinalBF';

A = dir([projfolder '\data\*.mat']);
for i = 1:length(A)
    ix = strfind(A(i).name,freq);
    station = A(i).name(1:ix-1);
    sel_data.results_dir = [projfolder '/' station '/'];
    sel_data.input_file = [projfolder '/data/' A(i).name];
    sel_data.modpath = [];
    sel_data.bfflag = 1;
    savepath = sel_data.results_dir;
    if exist([savepath '/statresults' addname '.mat'],'file')
        continue
    end
    plotallHDMCMCres2(sel_data);
end