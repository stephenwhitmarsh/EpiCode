
fname = 2711-p1-multifile-mHaT2_7

%% Add path
restoredefaultpath

if isunix
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/hspike/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/development/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/fieldtrip/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/DBSCAN/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/releaseDec2015/binaries/ 
end
if ispc
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\hspike
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\development
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\fieldtrip\
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\DBSCAN\
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\fieldtrip
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\MatlabImportExport_v6.0.0 % to read neuralynx files faster
end

ft_defaults

feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx

% get settings
config = hspike_setparams;

% get the artefacts
[MuseStruct_orig{ipatient}] = readMuseMarkers(config{ipatient}, false);

[marker{ipatient}, hypnogram{ipatient}] = hypnogramStats(config{ipatient}, MuseStruct_orig{ipatient}, false);

% trim files to only those within a hypnogram
MuseStruct_trimmed = MuseStruct_orig;
for ipart = 1 : 3
    sel = hypnogram{ipatient}.directory(hypnogram{ipatient}.part == ipart);
    first = find(strcmp(config{ipatient}.directorylist{ipart}, sel(1)));
    last = find(strcmp(config{ipatient}.directorylist{ipart}, sel(end)));
    config{ipatient}.directorylist{ipart} = config{ipatient}.directorylist{ipart}(first:last);
    MuseStruct_trimmed{ipatient}{ipart} = MuseStruct_trimmed{ipatient}{ipart}(first:last);
    if size(config{ipatient}.directorylist{ipart}, 2) > 7
        config{ipatient}.directorylist{ipart} = config{ipatient}.directorylist{ipart}(end-6:end);
        MuseStruct_trimmed{ipatient}{ipart} = MuseStruct_trimmed{ipatient}{ipart}(end-6:end);
    end
end

% write data concatinated for SC, artefacts, and output sampleinfo per file
writeSpykingCircus(config{ipatient}, MuseStruct_trimmed{ipatient}, true, true);

cfg = config{1};

st_FieldSelection(1) = 1; %timestamps
st_FieldSelection(2) = 1; %Channel Numbers
st_FieldSelection(3) = 1; %sample freq
st_FieldSelection(4) = 1; %Number of Valid Samples
st_FieldSelection(5) = 1; %samples
st_FieldSelection(6) = 1; %header (for creating .ncs only)
s_ExtractHeader      = 1;
s_ExtractMode        = 1; %1 if all samples
v_ModeArray          = []; %[] if all, 2 elements if range
s_AppendToFileFlag   = 0;
s_ExportMode         = 1; %1 if a
v_ExportModeVector   = [];

fname = fullfile(cfg.datasavedir, 'test_in.ncs');

[v_Timestamp, v_ChanNum, v_SampleFrequency, v_NumValSamples, v_Samples, st_Header] = ...
    Nlx2MatCSC_v3(fname, st_FieldSelection(1:5), s_ExtractHeader, s_ExtractMode, v_ModeArray);

st_FieldSelection(1) = 0; %timestamps


[v_ChanNum, v_SampleFrequency, v_NumValSamples, v_Samples, st_Header] = ...
    Nlx2MatCSC_v3(fname, st_FieldSelection(1:5), s_ExtractHeader, s_ExtractMode, v_ModeArray);


fname = fullfile(cfg.datasavedir, 'test_uit.ncs');

if isunix
    v_NumRecs = length(v_Timestamp);
    Mat2NlxCSC(fname, s_AppendToFileFlag, s_ExportMode, v_ExportModeVector, v_NumRecs, st_FieldSelection, ...
        v_Timestamp, v_ChanNum, v_SampleFrequency,  v_NumValSamples, v_Samples, st_Header);
end

[v_Timestamp, v_ChanNum, v_SampleFrequency, v_NumValSamples, v_Samples, st_Header] = ...
    Nlx2MatCSC_v3(fname, st_FieldSelection(1:5), s_ExtractHeader, s_ExtractMode, v_ModeArray);


















fin                         = '\\lexport\iss01.epimicro\patients\shared\raw_output\pat_02711_1193\eegcontinuous\02711_2019-04-18_11-04\02711_2019-04-18_11-04_mHaT2_2.ncs';
fout                        = '\\lexport\iss01.charpier\analyses\stephen.whitmarsh\data\test.ncs';

hdrin                       = ft_read_header(fin);
dat                         = ft_read_data(fin);

hdr                         = [];
hdr.Fs                      = hdrin.Fs;
hdr.nSamples                = size(dat, 2);
hdr.nSamplePre              = 0;
hdr.nChans                  = 1;
hdr.FirstTimeStamp          = 0;
hdr.TimeStampPerSample      = hdrin.TimeStampPerSample;
hdr.label{1}                = 'micro_label';
ft_write_data(fout, dat, 'chanindx', 1, 'dataformat', 'neuralynx_ncs', 'header', hdr);


hdrin2                      = ft_read_header(fout);
dat2                        = ft_read_data(fout);

fcon = '\\lexport\iss01.charpier\analyses\stephen.whitmarsh\data\hspike\2711\p1\2711-p1-multifile-1.ncs';
hdrcon                      = ft_read_header(fcon);






f{1} = '\\lexport\iss01.epimicro\patients\shared\raw_output\pat_02711_1193\eegcontinuous\02711_2019-04-18_11-04\02711_2019-04-18_11-04_mHaT2_2.ncs';
dat1 = ft_read_neuralynx_interp(f);
f{1} = '\\lexport\iss01.epimicro\patients\raw\pat_02711_1193\eeg\02711_2019-04-18_13-04\02711_2019-04-18_13-04_mHaT2_2.ncs';
dat2 = ft_read_neuralynx_interp(f);
dat = [dat1.trial{1} dat2.trial{1}];

clear f
f    = '\\lexport\iss01.epimicro\patients\shared\raw_output\pat_02711_1193\eegcontinuous\02711_2019-04-18_11-04\02711_2019-04-18_11-04_mHaT2_7.ncs';
dat1 = ft_read_data(f);
f    = '\\lexport\iss01.charpier\analyses\stephen.whitmarsh\datatest\02711_2019-04-18_10-29\02711_2019-04-18_10-29_mHaT2_7.ncs';
dat2 = ft_read_data(f);

% f    = '\\lexport\iss01.epimicro\patients\raw\pat_02711_1193\eeg\02711_2019-04-18_13-04\02711_2019-04-18_13-04_mHaT2_2.ncs';
% dat2 = ft_read_data(f);



dat = [dat1 dat2];
fout = '\\lexport\iss01.charpier\analyses\stephen.whitmarsh\data\test.ncs';
ft_write_data(fout, dat, 'chanindx', 1, 'dataformat', 'neuralynx_ncs', 'header', hdr);
u68
f{1} = fout;
datcon   = ft_read_neuralynx_interp(f);

fname = '\\lexport\iss01.charpier\analyses\stephen.whitmarsh\data\hspike\2711\p1\2711-p1-multifile-mHaT2_7.ncs';



fname = '/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/data/hspike/2711/p1/2711-p1-multifile-mHaT2_7.ncs';
hdr = ft_read_header(fname);


fname = '/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/data/hspike/2711/p2/2711-p2-multifile-mHaT2_7.ncs';
hdr = ft_read_header(fname);

fname = '/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/data/hspike/2711/p3/2711-p3-multifile-mHaT2_7.ncs';
hdr = ft_read_header(fname);




