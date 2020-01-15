addpath /network/lustre/iss01/home/laurent.bailly/fieldtrip/
ft_defaults
addpath /network/lustre/iss01/epimicro/patients/shared/Laurent/Scripts_matlab/

% addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip/

ft_defaults

%% Settings
vect_pat            = [2067,  2098,  2141,  2161,  2222,  2256,  2230, 2316,  2306,  2357,  2379,  2578,  2599,  2614,   2619, 2651  ];     % number of patients
diff_pat            = [false, false, false, false, false, false, true, false, false, false, false, false, false, false,  true, false ];     % Does patient have differnt patterns (in sleep)?
nrelectrodes(1,1,:) = [1,     1,     2,     2,     1,     2,     2,    3,     1,     1,     2,     1,     1,     4,      3,    1     ];     % Macro sleep/sommeil
nrelectrodes(1,2,:) = [1,     1,     2,     2,     1,     2,     2,    3,     1,     1,     2,     1,     1,     4,      3,    1     ];     % Macro wake/veille
nrelectrodes(2,1,:) = [0,     0,     1,     1,     1,     0,     1,    0,     1,     1,     2,     0,     1,     3,      2,    1     ];     % Micro sleep/sommeil
nrelectrodes(2,2,:) = [0,     0,     1,     1,     1,     0,     1,    0,     1,     0,     2,     0,     1,     3,      2,    1     ];     % Micro wake/veille

% location code: 1 = hippocampus, 2 = lesion, 3 = other, 0 = no electrode
location(1,1,:)     = [1,     1,     1,     1,     1,     3,     2,    1,     1,     1,     1,     1,     1,     1,      1,    1     ];     % location code Macro, electrode 1    
location(1,2,:)     = [0,     0,     1,     1,     0,     3,     1,    1,     0,     0,     1,     0,     0,     2,      1,    0     ];     % location code Macro, electrode 2    
location(1,3,:)     = [0,     0,     0,     0,     0,     0,     0,    3,     0,     0,     0,     0,     0,     2,      2,    0     ];     % location code Macro, electrode 3    
location(1,4,:)     = [0,     0,     0,     0,     0,     0,     0,    3,     0,     0,     0,     0,     0,     3,      0,    0     ];     % location code Macro, electrode 4    
location(2,1,:)     = [0,     1,     2,     2,     1,     2,     2,    3,     1,     1,     2,     0,     0,     4,      3,    0     ];     % location code Micro, electrode 1    
location(2,2,:)     = [0,     0,     1,     1,     1,     0,     2,    0,     1,     1,     1,     1,     1,     1,      1,    1     ];     % location code Micro, electrode 2    
location(2,3,:)     = [0,     0,     0,     0,     0,     0,     2,    0,     0,     0,     1,     0,     0,     2,      1,    0     ];     % location code Micro, electrode 3    
location(2,4,:)     = [0,     0,     0,     0,     0,     0,     2,    0,     0,     0,     0,     0,     0,     2,      0,    0     ];     % location code Micro, electrode 4    

electrode           = {'macro', 'micro'};        % electrode type
period              = {'sommeil', 'veille'};     % period
p_som               = {'P1', 'P2'};              % if differents periodes
datadir             = '/network/lustre/iss01/epimicro/patients/shared/Laurent/RIPPLELAB_Paris/Analysis';
analysisdatadir     = '/network/lustre/iss01/epimicro/patients/shared/Laurent/RIPPLELAB_Paris/Analysis/Results'; % directory to save file

% limits of time axis for all plotting (and TFR analysis);
xlim_TFR = [-0.1,0.1];
xlim_dat = [-0.4,0.4];

% write Fieldtrip structures from RippleLab data - skip if done already
% ripples_getdata

%% load data

fout = fullfile(analysisdatadir,'eventdat');
load(fout,'eventdat')

hasdat = ind2sub(size(eventdat),~cellfun('isempty', eventdat));

% combine data to average micro-macro
cfg                         = [];
cfg.keepsampleinfo          = 'no';

ER_ripple_multi_macro       = ft_appenddata(cfg,eventdat{hasdat(:,1) & electrodetype' == 1,1});
ER_ripple_single_macro      = ft_appenddata(cfg,eventdat{hasdat(:,2) & electrodetype' == 1,2});
ER_fastripple_multi_macro   = ft_appenddata(cfg,eventdat{hasdat(:,4) & electrodetype' == 1,4});
ER_fastripple_single_macro  = ft_appenddata(cfg,eventdat{hasdat(:,3) & electrodetype' == 1,3});
ER_ripple_multi_micro       = ft_appenddata(cfg,eventdat{hasdat(:,1) & electrodetype' == 2,1});
ER_ripple_single_micro      = ft_appenddata(cfg,eventdat{hasdat(:,2) & electrodetype' == 2,2});
ER_fastripple_multi_micro   = ft_appenddata(cfg,eventdat{hasdat(:,4) & electrodetype' == 2,4});
ER_fastripple_single_micro  = ft_appenddata(cfg,eventdat{hasdat(:,3) & electrodetype' == 2,3});

% loop over locations
for iloc = 1 : 3 
    try ER_ripple_multi_macro_loc{iloc}         = ft_appenddata(cfg,eventdat{hasdat(:,1) & electrodetype' == 1 & eleclocation' == iloc,1}); catch end
    try ER_ripple_single_macro_loc{iloc}        = ft_appenddata(cfg,eventdat{hasdat(:,2) & electrodetype' == 1 & eleclocation' == iloc,2}); catch end
    try ER_fastripple_multi_macro_loc{iloc}     = ft_appenddata(cfg,eventdat{hasdat(:,4) & electrodetype' == 1 & eleclocation' == iloc,4}); catch end
    try ER_fastripple_single_macro_loc{iloc}    = ft_appenddata(cfg,eventdat{hasdat(:,3) & electrodetype' == 1 & eleclocation' == iloc,3}); catch end
    try ER_ripple_multi_micro_loc{iloc}         = ft_appenddata(cfg,eventdat{hasdat(:,1) & electrodetype' == 2 & eleclocation' == iloc,1}); catch end
    try ER_ripple_single_micro_loc{iloc}        = ft_appenddata(cfg,eventdat{hasdat(:,2) & electrodetype' == 2 & eleclocation' == iloc,2}); catch end
    try ER_fastripple_multi_micro_loc{iloc}     = ft_appenddata(cfg,eventdat{hasdat(:,4) & electrodetype' == 2 & eleclocation' == iloc,4}); catch end
    try ER_fastripple_single_micro_loc{iloc}    = ft_appenddata(cfg,eventdat{hasdat(:,3) & electrodetype' == 2 & eleclocation' == iloc,3}); catch end
end


% load data of examples
cfg             = [];
cfg.resamplefs  = 1000;
temp            = load(fullfile(analysisdatadir,'2098-sommeil-macro-1-data.mat')); cfg.trials = 19;  ex1 = ft_resampledata(cfg,temp.dat);
temp            = load(fullfile(analysisdatadir,'2222-sommeil-macro-1-data.mat')); cfg.trials = 101; ex2 = ft_resampledata(cfg,temp.dat);
temp            = load(fullfile(analysisdatadir,'2098-sommeil-macro-1-data.mat')); cfg.trials = 7;   ex3 = ft_resampledata(cfg,temp.dat);
temp            = load(fullfile(analysisdatadir,'2098-sommeil-macro-1-data.mat')); cfg.trials = 43;  ex4 = ft_resampledata(cfg,temp.dat);

% TFR analysis (often problem with java)
cfg             = [];
cfg.output      = 'pow';
cfg.channel     = 'all';
cfg.method      = 'mtmconvol';
cfg.taper       = 'hanning';
cfg.pad         = 'nextpow2';
cfg.foi         = 1:5:500;                         % try with steps of 1 hertz:
cfg.t_ftimwin   = 7./cfg.foi;                       % 7 cycles per time window, can try with 5
cfg.toi         = xlim_TFR(1):0.001:xlim_TFR(2);    % can try changing xlim_TFR for more/less zoom, and the stepsize for more smooth estimates over time (but 0.001 is already pretty small)

% macro electrodes
TFR_ripple_multi_macro       = ft_freqanalysis(cfg, ER_ripple_multi_macro);
TFR_ripple_single_macro      = ft_freqanalysis(cfg, ER_ripple_single_macro);
TFR_fastripple_multi_macro   = ft_freqanalysis(cfg, ER_fastripple_multi_macro);
TFR_fastripple_single_macro  = ft_freqanalysis(cfg, ER_fastripple_single_macro);

% micro electrodes
TFR_ripple_multi_micro       = ft_freqanalysis(cfg, ER_ripple_multi_micro);
TFR_ripple_single_micro      = ft_freqanalysis(cfg, ER_ripple_single_micro);
TFR_fastripple_multi_micro   = ft_freqanalysis(cfg, ER_fastripple_multi_micro);
TFR_fastripple_single_micro  = ft_freqanalysis(cfg, ER_fastripple_single_micro);

% loop over locations
for iloc = 1 : 3
    % macro electrodes
    try TFR_ripple_multi_macro_loc{iloc}       = ft_freqanalysis(cfg, ER_ripple_multi_macro_loc{iloc}); catch end
    try TFR_ripple_single_macro_loc{iloc}      = ft_freqanalysis(cfg, ER_ripple_single_macro_loc{iloc}); catch end
    try TFR_fastripple_multi_macro_loc{iloc}   = ft_freqanalysis(cfg, ER_fastripple_multi_macro_loc{iloc}); catch end
    try TFR_fastripple_single_macro_loc{iloc}  = ft_freqanalysis(cfg, ER_fastripple_single_macro_loc{iloc}); catch end
    
    % micro electrodes
    try TFR_ripple_multi_micro_loc{iloc}       = ft_freqanalysis(cfg, ER_ripple_multi_micro_loc{iloc}); catch end
    try TFR_ripple_single_micro_loc{iloc}      = ft_freqanalysis(cfg, ER_ripple_single_micro_loc{iloc}); catch end
    try TFR_fastripple_multi_micro_loc{iloc}   = ft_freqanalysis(cfg, ER_fastripple_multi_micro_loc{iloc}); catch end
    try TFR_fastripple_single_micro_loc{iloc}  = ft_freqanalysis(cfg, ER_fastripple_single_micro_loc{iloc}); catch end
end


% examples
ex1TFR          = ft_freqanalysis(cfg, ex1);
ex2TFR          = ft_freqanalysis(cfg, ex2);
ex3TFR          = ft_freqanalysis(cfg, ex3);
ex4TFR          = ft_freqanalysis(cfg, ex4);

% save TFR
fout = fullfile(analysisdatadir,'TFR');
save(fout,'TFR*','ex*','-v7.3');

% load TFR
% fout = fullfile(analysisdatadir,'TFR');
% load(fout,'TFR*','ex*');


%%%%%%%%%%%%%%%%
%%% PLOTTING %%%
%%%%%%%%%%%%%%%%

%% Examples

fig                 = figure;
colormap jet % instead of default: parula

cfg                 = [];
cfg.xlim            = xlim_dat;

subplot(2,4,1); ft_singleplotER(cfg,ex1); title('Ripple');
subplot(2,4,2); ft_singleplotER(cfg,ex2); title('Ripple Multi');
subplot(2,4,3); ft_singleplotER(cfg,ex3); title('FastRipple');
subplot(2,4,4); ft_singleplotER(cfg,ex4); title('Fast-on-Ripple');

cfg                 = [];
cfg.baseline        = 'yes';
cfg.baselinetype    = 'relative';
cfg.xlim            = xlim_TFR;
cfg.zlim            = 'maxabs';
cfg.colorbar        = 'no';

% TFR average
subplot(2,4,5); ft_singleplotTFR(cfg,ex1TFR); colorbar('southoutside'); title('');
subplot(2,4,6); ft_singleplotTFR(cfg,ex2TFR); colorbar('southoutside'); title(''); 
subplot(2,4,7); ft_singleplotTFR(cfg,ex3TFR); colorbar('southoutside'); title(''); 
subplot(2,4,8); ft_singleplotTFR(cfg,ex4TFR); colorbar('southoutside'); title(''); 

% print to file
fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', fullfile(analysisdatadir,'examples.pdf'),'-r300');


%% Plot average micro

fig = figure;
colormap jet % instead of default: parula

cfg = [];
cfg.xlim = xlim_dat; % change it at the top of the script!

subplot(2,4,1); ft_singleplotER(cfg,ER_ripple_single_micro);       title('Ripple');
subplot(2,4,2); ft_singleplotER(cfg,ER_ripple_multi_micro);        title('Ripple Multi');
subplot(2,4,3); ft_singleplotER(cfg,ER_fastripple_single_micro);   title('FastRipple');
subplot(2,4,4); ft_singleplotER(cfg,ER_fastripple_multi_micro);    title('Fast-on-Ripple');

cfg = [];
cfg.baseline        = 'yes';
cfg.baselinetype    = 'relative';
cfg.xlim            = xlim_TFR;
cfg.zlim            = 'maxabs';
cfg.colorbar        = 'no';
cfg.title           = '';

% TFR average
subplot(2,4,5); ft_singleplotTFR(cfg,TFR_ripple_single_micro);     colorbar('southoutside');  title('');
subplot(2,4,6); ft_singleplotTFR(cfg,TFR_ripple_multi_micro);      colorbar('southoutside');  title('');
subplot(2,4,7); ft_singleplotTFR(cfg,TFR_fastripple_single_micro); colorbar('southoutside');  title('');
subplot(2,4,8); ft_singleplotTFR(cfg,TFR_fastripple_multi_micro);  colorbar('southoutside');  title('');

% print to file
fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', fullfile(analysisdatadir,'averages_micro.pdf'),'-r300');


%% Plot average macro

fig = figure;
%colormap jet % instead of default: parula

cfg = [];
cfg.xlim = xlim_dat; % change it at the top of the script!

subplot(2,4,1); ft_singleplotER(cfg,ER_ripple_single_macro);       title('Ripple');
subplot(2,4,2); ft_singleplotER(cfg,ER_ripple_multi_macro);        title('Ripple Multi');
subplot(2,4,3); ft_singleplotER(cfg,ER_fastripple_single_macro);   title('FastRipple');
subplot(2,4,4); ft_singleplotER(cfg,ER_fastripple_multi_macro);    title('Fast-on-Ripple');

cfg = [];
cfg.baseline        = 'yes';
cfg.baselinetype    = 'relative';
cfg.xlim            = xlim_TFR; % change it at the top of the script!
cfg.zlim            = 'maxabs';
cfg.colorbar        = 'no';
cfg.title           = '';

% TFR average
subplot(2,4,5); ft_singleplotTFR(cfg,TFR_ripple_single_macro);     colorbar('southoutside'); title(''); 
subplot(2,4,6); ft_singleplotTFR(cfg,TFR_ripple_multi_macro);      colorbar('southoutside'); title(''); 
subplot(2,4,7); ft_singleplotTFR(cfg,TFR_fastripple_single_macro); colorbar('southoutside'); title(''); 
subplot(2,4,8); ft_singleplotTFR(cfg,TFR_fastripple_multi_macro);  colorbar('southoutside'); title(''); 

% print to file
fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', fullfile(analysisdatadir,'averages_macro.pdf'),'-r300');




%% Plot average macro for different locations

for iloc = 1 : 3
    
    fig = figure;
    %colormap jet % instead of default: parula
    
    cfg                 = [];
    cfg.xlim            = xlim_TFR; % change it at the top of the script!
    
    subplot(2,4,1); try ft_singleplotER(cfg,ER_ripple_single_macro_loc{iloc});       title('Ripple');           catch end
    subplot(2,4,2); try ft_singleplotER(cfg,ER_ripple_multi_macro_loc{iloc});        title('Ripple Multi');     catch end
    subplot(2,4,3); try ft_singleplotER(cfg,ER_fastripple_single_macro_loc{iloc});   title('FastRipple');       catch end
    subplot(2,4,4); try ft_singleplotER(cfg,ER_fastripple_multi_macro_loc{iloc});    title('Fast-on-Ripple');   catch end
    
    cfg = [];
    cfg.baseline        = 'yes';
    cfg.baselinetype    = 'relative';
    cfg.xlim            = xlim_dat; % change it at the top of the script!
    cfg.zlim            = 'maxabs';
    cfg.colorbar        = 'no';
    cfg.title           = '';
    
    % TFR average
    subplot(2,4,5); try ft_singleplotTFR(cfg,TFR_ripple_single_macro_loc{iloc});     colorbar('southoutside'); title(''); catch end
    subplot(2,4,6); try ft_singleplotTFR(cfg,TFR_ripple_multi_macro_loc{iloc});      colorbar('southoutside'); title(''); catch end
    subplot(2,4,7); try ft_singleplotTFR(cfg,TFR_fastripple_single_macro_loc{iloc}); colorbar('southoutside'); title(''); catch end
    subplot(2,4,8); try ft_singleplotTFR(cfg,TFR_fastripple_multi_macro_loc{iloc});  colorbar('southoutside'); title(''); catch end
    
    % print to file
    fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(analysisdatadir,['averages_macro_loc',num2str(iloc),'.pdf']),'-r300');
end



%% Plot average micro for different locations

for iloc = 1 : 3
    
    fig = figure;
    %colormap jet % instead of default: parula
    
    cfg                 = [];
    cfg.xlim            = xlim_TFR; % change it at the top of the script!
    
    subplot(2,4,1); try ft_singleplotER(cfg,ER_ripple_single_micro_loc{iloc});       title('Ripple');           catch end
    subplot(2,4,2); try ft_singleplotER(cfg,ER_ripple_multi_micro_loc{iloc});        title('Ripple Multi');     catch end
    subplot(2,4,3); try ft_singleplotER(cfg,ER_fastripple_single_micro_loc{iloc});   title('FastRipple');       catch end
    subplot(2,4,4); try ft_singleplotER(cfg,ER_fastripple_multi_micro_loc{iloc});    title('Fast-on-Ripple');   catch end
    
    cfg = [];
    cfg.baseline        = 'yes';
    cfg.baselinetype    = 'relative';
    cfg.xlim            = xlim_dat; % change it at the top of the script!
    cfg.zlim            = 'maxabs';
    cfg.colorbar        = 'no';
    cfg.title           = '';
    
    % TFR average
    subplot(2,4,5); try ft_singleplotTFR(cfg,TFR_ripple_single_micro_loc{iloc});     colorbar('southoutside'); title(''); catch end
    subplot(2,4,6); try ft_singleplotTFR(cfg,TFR_ripple_multi_micro_loc{iloc});      colorbar('southoutside'); title(''); catch end
    subplot(2,4,7); try ft_singleplotTFR(cfg,TFR_fastripple_single_micro_loc{iloc}); colorbar('southoutside'); title(''); catch end
    subplot(2,4,8); try ft_singleplotTFR(cfg,TFR_fastripple_multi_micro_loc{iloc});  colorbar('southoutside'); title(''); catch end
    
    % print to file
    fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(analysisdatadir,['averages_micro_loc',num2str(iloc),'.pdf']),'-r300');
end


