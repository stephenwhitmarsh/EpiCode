try %en local
    scriptpath = matlab.desktop.editor.getActiveFilename;
catch %cluster
    scriptpath = mfilename('fullpath');
end

epicodepath = [fileparts(fileparts(fileparts(scriptpath))), filesep];

addpath (genpath([epicodepath,'development']))
addpath (genpath([epicodepath,'shared']))
addpath (genpath([epicodepath,'external']))
addpath (genpath([epicodepath,'templates']))
addpath (genpath([epicodepath, 'wod']))
addpath (genpath([epicodepath, 'dtx']))
addpath (genpath([epicodepath, 'wod',filesep,'wod_functions']))
addpath (genpath([epicodepath, 'wod',filesep,'Intra']))

if ispc
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
    addpath \\lexport\iss01.charpier\echanges\scripts-paul\Spike2_vers_MATLAB
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
    addpath /network/lustre/iss01/charpier/echanges/scripts-paul/Spike2_vers_MATLAB

end

ft_defaults

config = Intra_setparams;


chanlist= ["Vm","EEG-S1-L"];%,"Im"]; %chanlist Intra

for iprot= 13:16%size(config,2)
    datapath = char(fullfile(config{iprot}.rawdir,config{iprot}.directorylist{1}));
    matlabpath= fullfile(fileparts(datapath),'matlab_structures');
    [folder, file, extension] = fileparts(datapath);

    CEDStruct = readCEDevents(datapath);
    
    for channame = chanlist
        data = readCEDwaveforms(datapath, channame);
        
        

        cfgtemp= [];
        cfgtemp.resamplefs      = config{iprot}.Intra.resamplefs;
        cfgtemp.detrend         = 'no';
        data                    = ft_resampledata(cfgtemp,data);
        

        if ~isfolder(matlabpath)
            mkdir(matlabpath);
        end

        fname=fullfile(matlabpath,sprintf('%s_%s.mat',config{iprot}.directorylist{1}{1},channame));
        save(fname,'data');
        clear fname
        fname=fullfile(matlabpath,sprintf('%s_events.mat',config{iprot}.directorylist{1}{1}));
        save(fname,'CEDStruct');
    end %channame
end %iprot
