%% Setting parameters

%
% cfg.patientdir            = [path to patient data, one above /eeg and /eegmicromed]
% cfg.hyp.micromedchannel   = [micromed channel to use for alignment]
% cfg.imagesavedir          = [path to directory in which to save plots]
% cfg.prefix                = [string to prefix output figure, e.g. patient ID]
% cfg.hyp.backupdir         = [directory to backup Muse marker files]
%
% Example:
%
% cfg                       = [];
% cfg.patientdir            = '/network/lustre/iss01/epimicro/patients/raw/pat_02711_1193';
% cfg.hyp.micromedchannel   = 'F3p6';
% cfg.imagesavedir          = '/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/images/hspike';
% cfg.prefix                = 'P1-';
% cfg.hyp.backupdir         = '/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/markerbackup';
% export_hypnogram(cfg);


function [cfg] = sleep_setparams(cfg, patID)

if ~exist(cfg.patientdir, 'dir')
    mkdir(cfg.patientdir)
end
if ~exist(cfg.imagesavedir, 'dir')
    mkdir(cfg.imagesavedir)
end
if ~exist(cfg.hyp.backupdir, 'dir')
    mkdir(cfg.hyp.backupdir)
end

switch patID
    case 'pat_02614_1073'
        cfg.hyp.micromedchannel   = 'Caid1';
        cfg.prefix                = [patID '_'];
    case 'pat_02619_1078'
        cfg.hyp.micromedchannel   = 'AmT21';
        cfg.prefix                = [patID '_'];
    case 'pat_02651_1127'
        cfg.hyp.micromedchannel   = 'ACal2';
        cfg.prefix                = [patID '_'];
    case 'pat_02660_1136'
        cfg.hyp.micromedchannel   = 'Am2g2';
        cfg.prefix                = [patID '_'];
    case 'pat_02680_1158'
        cfg.hyp.micromedchannel   = 'AmT22';
        cfg.prefix                = [patID '_'];
    case 'pat_02689_1168'
        cfg.hyp.micromedchannel   = 'Car2';
        cfg.prefix                = [patID '_'];
    case 'pat_02711_1193'
        cfg.hyp.micromedchannel   = 'F3p6';
        cfg.prefix                = [patID '_'];
    case 'pat_02718_1201'
        cfg.hyp.micromedchannel   = 'HaT12';
        cfg.prefix                = [patID '_'];
    case 'pat_02757_1244'
        cfg.hyp.micromedchannel   = 'F1a2';
        cfg.prefix                = [patID '_'];
        
end
