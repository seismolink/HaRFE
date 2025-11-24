function start_harfe(func,varargin)

% Initialization of HaRFE
% Probabilistic Inversion of Receiver Functions for Anisotropic Structure 
% based on Azimuthaly Harmonic Variations.

% input
%   path2opts: path to file with options for running RFprep
%   func:
%       +version:  	        print version of current software
%       +check:             check if a newer version is available
%       +prepare:           prepare project for inversion

% Copyright 2025 F.Link

close all
clc

if nargin < 1
    func = '';
end

% clear all
clearvars -except path2opts func varargin
clear functions
clear global

curr_version = '1.0.0';
program_url = 'https://github.com/seismolink/HaRFE/';

warning('off','all')

txt = strvcat({ ...
    'HaRFE - Copyright 2025 F. Link', ...
    'This program is licensed under the Apache License, Version 2.0 ', ...
    '(the "License"); you may not use this file except in compliance ', ...
    'with the License. You may obtain a copy of the License at:', ...
    ' ', ...
    '    http://www.apache.org/licenses/LICENSE-2.0', ...
    ' ', ...
    'Unless required by applicable law or agreed to in writing, software ', ...
    'distributed under the License is distributed on an "AS IS" BASIS, ', ...
    'WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or ', ...
    'implied. See the License for the specific language governing ', ...
    'permissions and limitations under the License.', ...
    ' ', ...
    'If the use or modification of this program contributes to a ', ...
    'scientific publication, please consider citing the accompanying ', ...
    'publications.' ...
});

disp(txt);

% get current folder name
cur_dir = pwd;
[~, df, ~] = fileparts(cur_dir); 

% add all sub directories to search path
addpath(genpath(['../',df]));

% Check which actions are requested
vrs_flag = 0;
chck_flag = 0;
prep_flag = 0;
if strcmp(func,'+version')
    vrs_flag = 1;
end
if strcmp(func,'+check')
    chck_flag = 1;
end
if strcmp(func,'+prepare')
    prep_flag = 1;
end

if ~vrs_flag && ~chck_flag && ~prep_flag
    disp(' ')
    disp(' ')
    disp('HaRFE is a program for the Probabilistic Inversion of Receiver Functions for Anisotropic Structure.')
    disp('---------------------------------------------------------------------------------------------------')
    disp('input')
    disp('  func:')
    disp('      +version:  	    print version of current software')
    disp('      +check:         check if a newer version is available')
    disp('      +prepare:       prepare project for inversion')
end


% perform requested actions
if vrs_flag
    disp(['Current version of HaRFE: ' curr_version]);
    return
end
if chck_flag
    disp(['Current version of HaRFE: ' curr_version]);
    disp('Checking online for recent version. Please wait...')
    if ispc
        request_url = 'wget64';
    else
        request_url = 'wget';
    end
    request_url = strcat(request_url,' "',program_url,'releases/latest" -O temp.txt');
    [~,~] = system(request_url);
    fid = fopen('temp.txt','rt');
    line = fgetl(fid);
    while ischar(line)
        if contains(line,'<title>')
            disp(line(10:end-8));
            if ~contains(line,curr_version)
                fprintf(strcat('It is likely that a new version is available.\n',...
                    'Visit %s to download the newest version.\n'),program_url);
            else
                disp('Current version is the most recent.')
            end
            break;
        end
        line = fgetl(fid);
    end
    fclose(fid);
    delete('temp.txt')
    return
end
if prep_flag
    prepareINV;
end
end