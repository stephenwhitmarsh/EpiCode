function [SpikeIndex, ChanId, SpikeFV] = mDetectSpike(varargin)
% This function detects spikes in the EEG. The input parameters to the function are
% (1) DATA = EEG signal of M x N dimension where M = data length and N = # of
% channels. 
% (2) Fs = samping rate
% (3) BlockSize = length of data block to be processed in minutes
% (default = 1 minute block).
% (4) SCALE = scale the EEG amplitude in percentage. Default is 70%.
% (5) STDCoeff = distance from the mean (measure of spread across the mean
% default = 4 times the standard deviation from the mean).
% (6) DetThresholds = specify the detection thresholds which contains
% thresholds for the left half-wave slope, right half-wave slope, total
% amplitude of the spike (left half-wave + right half-wave), duration of
% left half-wave and right half-wave in ms. Default values for the  
% DetThresholds = [7; 7; 600; 10 ; 10];
% (7) FilterSpec contains the filter specification. The cut-off frequency
% for the first high-pass and low-pass filter, followed by cut-off
% frequency for the second high-pass and low-pass filter. Filter order is
% kept fixed.
% (8) TroughSearch is distance in ms to search for a trough on each side of
% a detected peak. Default width is 40 ms.
% The function returns returns time-instance of detected spike events (SpikeIndex),
% detection channel (ChanId), and spike feature (SpikeFV). The SpikeFV is of 
% dimension P x 8, where P represents numbers of spikes detected in the data: 
% (1)index of peak in block, 
% (2)index of peak in whole dataset, 
% (3)Value at peak of spike, 
% (4)Left half-wave amplitude,
% (5)Left half-wave duration, 
% (6)Left half-wave slope, 
% (7)Right half-wave amplitude, 
% (8)Right half-wave duration, and
% (9)Right half-wave slope.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Usage: [SpikeIndex, ChanId, SpikeFV] = mDetectSpike(Data, Fs, BlockSize,
% SCALE, STDCoeff, DetThresholds, FilterSpec, TroughSearch);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CodeID:           SPKDT v1.0.4
% Author:           Danielle T. Barkmeier
% Modified by:      Rajeev Yadav
% Email:            rajeevyadav@gmail.com
% Dated:            September, 11th 2011
% Rev. Update:
%                   1. Fixed bug in filter specification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Screen user input parameters
if nargin < 2
    disp('Insufficient number of argumens');
    return;
else switch nargin
        case 2
            Data = varargin{1};
            Fs = varargin{2};
            LS = 7;         % Left half-wave slope
            RS = 7;         % Right half-wave slope
            TAMP = 600;     % Totoal amplitude
            LD = 10;        % Left half-wave duration
            RD = 10;        % Right half-wave duration
            STDCoeff = 4;   % Chebyshev inequality coefficient (distance from centre point or mean)
            SCALE = 70;     % Scaling parameter
            BlockSize = 1;  % Data processing block size in minutes
            TroughSearch = 40; % distance in ms to search for a trough on each side of a detected peak
            DetThresholds = [ LS; RS; TAMP; LD; RD;];
            FilterSpec =  [20; 50; 1; 35;];
        case 3
            Data = varargin{1};
            Fs = varargin{2};
            BlockSize = varargin{3};
            LS = 7;         % Left half-wave slope
            RS = 7;         % Right half-wave slope
            TAMP = 600;     % Totoal amplitude
            LD = 10;        % Left half-wave duration
            RD = 10;        % Right half-wave duration
            STDCoeff = 4;   % Chebyshev inequality coefficient (distance from centre point or mean)
            SCALE = 70;     % Scaling parameter
            TroughSearch = 40;
            DetThresholds = [ LS; RS; TAMP; LD; RD;];
            FilterSpec =  [20; 50; 1; 35;];
        case 4
            Data = varargin{1};
            Fs = varargin{2};
            BlockSize = varargin{3};
            SCALE = varargin{4}; 
            LS = 7;         % Left half-wave slope
            RS = 7;         % Right half-wave slope
            TAMP = 600;     % Totoal amplitude
            LD = 10;        % Left half-wave duration
            RD = 10;        % Right half-wave duration
            STDCoeff = 4;   % Chebyshev inequality coefficient (distance from centre point or mean)
            TroughSearch = 40;
            FilterSpec =  [20; 50; 1; 35;];
            DetThresholds = [ LS; RS; TAMP; LD; RD;];
        case 5
            Data = varargin{1};
            Fs = varargin{2};
            BlockSize = varargin{3};
            SCALE = varargin{4};
            STDCoeff = varargin{5};
            LS = 7;         % Left half-wave slope
            RS = 7;         % Right half-wave slope
            TAMP = 600;     % Totoal amplitude
            LD = 10;        % Left half-wave duration
            RD = 10;        % Right half-wave duration
            TroughSearch = 40;
            DetThresholds = [ LS; RS; TAMP; LD; RD;];
            FilterSpec =  [20; 50; 1; 35;];
        case 6
            Data = varargin{1};
            Fs = varargin{2};
            BlockSize = varargin{3};
            SCALE = varargin{4};
            STDCoeff = varargin{5};
            DetThresholds = varargin{6};
            if length(DetThresholds) < 2
                LS = 7;         % Left half-wave slope
                RS = 7;         % Right half-wave slope
                TAMP = 600;     % Totoal amplitude
                LD = 10;        % Left half-wave duration
                RD = 10;        % Right half-wave duration
            elseif length(DetThresholds) == 5
                LS = DetThresholds(1); 
                RS = DetThresholds(2); 
                TAMP = DetThresholds(3); 
                LD = DetThresholds(4); 
                RD = DetThresholds(5);
            end
            TroughSearch = 40;
            FilterSpec =  [20; 50; 1; 35;];
        case 7
            Data = varargin{1};
            Fs = varargin{2};
            BlockSize = varargin{3};
            SCALE = varargin{4};
            STDCoeff = varargin{5};
            DetThresholds = varargin{6};
            if length(DetThresholds) < 2
                LS = 7;         % Left half-wave slope
                RS = 7;         % Right half-wave slope
                TAMP = 600;     % Totoal amplitude
                LD = 10;        % Left half-wave duration
                RD = 10;        % Right half-wave duration
            elseif length(DetThresholds) == 5
                LS = DetThresholds(1); 
                RS = DetThresholds(2); 
                TAMP = DetThresholds(3); 
                LD = DetThresholds(4); 
                RD = DetThresholds(5);
            end
            TroughSearch = 40;
            FilterSpec = varargin{7};
            if length(FilterSpec) < 2
                disp('using default filter specification');
                FilterSpec =  [20; 50; 1; 35;];
            end
        case 8
            Data = varargin{1};
            Fs = varargin{2};
            BlockSize = varargin{3};
            SCALE = varargin{4};
            STDCoeff = varargin{5};
            DetThresholds = varargin{6};
            if length(DetThresholds) < 2
                LS = 7;         % Left half-wave slope
                RS = 7;         % Right half-wave slope
                TAMP = 600;     % Totoal amplitude
                LD = 10;        % Left half-wave duration
                RD = 10;        % Right half-wave duration
            elseif length(DetThresholds) == 5
                LS = DetThresholds(1); 
                RS = DetThresholds(2); 
                TAMP = DetThresholds(3); 
                LD = DetThresholds(4); 
                RD = DetThresholds(5);
            end
            
            FilterSpec = varargin{7};
            if length(FilterSpec) < 2
                disp('using default filter specification');
                FilterSpec =  [20; 50; 1; 35;];
            end
                
            TroughSearch = varargin{8};
        otherwise
            error('Unexpected inputs');
    end
end

%% Check if input is empty or not. If empty, set defaults.
if ~exist('BlockSize','var')
    BlockSize = 1;
end
if ~exist('SCALE','var')
    SCALE = 70;
end
if ~exist('STDCoeff','var')
    STDCoeff = 4;
end
if ~exist('DetThresholds','var')
    LS = 7;         % Left half-wave slope
    RS = 7;         % Right half-wave slope
    TAMP = 600;     % Totoal amplitude
    LD = 10;        % Left half-wave duration
    RD = 10;        % Right half-wave duration
    DetThresholds = [ LS; RS; TAMP; LD; RD;];
end

if ~exist('FilterSpec','var')
    FilterSpec = [20; 50; 1; 35];
end
if length(FilterSpec) < 2
    disp('using default filter specification');
    FilterSpec =  [20; 50; 1; 35;];
end 

%% Filter specification
[bh1, ah1] = butter(2, FilterSpec(1)/(Fs/2), 'high');
[bl1, al1] = butter(4, FilterSpec(2)/(Fs/2), 'low');
[bh2, ah2] = butter(2, FilterSpec(3)/(Fs/2),  'high');
[bl2, al2] = butter(4, FilterSpec(4)/(Fs/2), 'low');
    
if isempty(TroughSearch)
    TroughSearch = 40;
end

%% Examine data and initialize parameters for processing
[M, N] = size(Data);            % where M = data size, N = number of channels
NumSecs = floor(M/Fs);                 % Get file length in seconds % added floor
ChanGroups = 1:N;               % Channel group
% Initialization
AllSpikes = [];
SpikeIndex = [];
ChanId = [];
SpikeFV  = [];

% Requires blocksize as input parameter
% round changes to floor by stephen
Blocks = floor(NumSecs/(BlockSize*60));

CurrentBlock = 1;               % Start index
for CurrentBlock = 1:Blocks
    if CurrentBlock == Blocks
        EEG = Data((CurrentBlock-1)*60*Fs+1:M,:); % get all data at the last data pointer
    else
        Idx = ((CurrentBlock-1)*60*Fs+1:(((CurrentBlock-1)*60*Fs+1)+ BlockSize*60*Fs));
        Idx = Idx(:);
        EEG = Data(Idx,:); 
    end
    EEG = single(EEG);          % convert data to single-precision
    c = 1/(1000/round(Fs));     % convert ms to # of data points 
    %% Detect artifact channel
    % Get average slope of each channel
    d1s = mean(abs(diff(EEG)))';  
    tmp = d1s;
    % Screen out any values more than twice the median slope
    tmp(find(tmp > median(tmp)*2)) = '';   
    
    % Detect artifact channels by those whose average slope are greater
    % than 10 standard deviations outside the mean of this screened data
    artifactChans = find(d1s > 10*std(tmp) + mean(tmp));
    
    % Remove them from consideration
    BlockChanNums = [];
    for kb = 1:size(ChanGroups,1)
        BlockChanNums = [BlockChanNums, setdiff(ChanGroups,artifactChans)];
    end


    %% Filter data in each channel for each block
    fEEG = [];
    EEG(:,BlockChanNums) = -EEG(:,BlockChanNums);% flip data
    fEEG(:,BlockChanNums) = filter(bh1, ah1, EEG(:,BlockChanNums));
    fEEG(:,BlockChanNums) = filter(bl1, al1, fEEG(:,BlockChanNums));
    % Filter original data
    EEG(:,BlockChanNums) = filter(bh2, ah2, EEG(:,BlockChanNums));
    EEG(:,BlockChanNums) = filter(bl2, al2, EEG(:,BlockChanNums));

    
    %% Scale channels
    % Default scaling = 70
    scale = SCALE/(median(mean(abs(EEG(:,ChanGroups)))));
    EEG(:,ChanGroups) = EEG(:,ChanGroups)*scale;
    
    %% Detect spikes
    % default STDCoeff = 4
   
    for Chan = BlockChanNums
        thresh = -mean(abs(fEEG(:,Chan))) - STDCoeff*std(abs(fEEG(:,Chan)));
        peaks = find(fEEG(100*c+1:end-100,Chan)<thresh);
        peaks = peaks+100*c+1;
        chanPeaks = zeros(1,9);
        for ki = 1:size(peaks,1)
            
            % added indexes by Stephen to fix when it goes over 
            starti = peaks(ki,1)-25*c;
            endi = peaks(ki,1)+25*c;
            if endi > size(fEEG,1)
                endi = size(fEEG,1);
            end
            [newPeakV newPeakI] = min(fEEG(starti:endi, Chan));
%             [newPeakV newPeakI] = min(fEEG(peaks(ki,1)-25*c:peaks(ki,1)+25*c, Chan));
            newPeakI = newPeakI + peaks(ki,1)-25*c-1;
            % See if this peak was already added to our list
            if(chanPeaks(end,1) ~= newPeakI)
                [spikeV spikeI] = max(EEG(newPeakI-20*c:newPeakI,Chan)); %Find peak
                spikeI = spikeI + newPeakI-20*c-1; 
                [leftV leftI] = min(EEG(spikeI-TroughSearch*c:spikeI,Chan));
                leftI = leftI + spikeI-TroughSearch*c-1; 
                
                % just to not break when looking beyond data in time
                try
                [rightV rightI] = min(EEG(spikeI:spikeI+TroughSearch*c,Chan));
                rightI = rightI + spikeI-1;
                catch 
                    
                end
                % Get amp, dur and slope for the left halfwave
                Lamp = spikeV-leftV;
                Ldur = (spikeI-leftI)/c;
                Lslope = Lamp/Ldur;
        
                % Get amp, dur and slope for the right halfwave
                Ramp = spikeV-rightV;
                Rdur = (rightI-spikeI)/c;
                Rslope = Ramp/Rdur;
                if(Lslope>LS&&Rslope>RS&&Lamp+Ramp>TAMP&&Ldur>LD&&Rdur>RD&& abs(chanPeaks(end,2)-spikeI)> 5*c)
                    chanPeaks = [chanPeaks; newPeakI,spikeI,spikeV,Lamp,Ldur,Lslope,Ramp,Rdur,Rslope];
                end
            end
        end
        % Next step ...
        chanPeaks = double(chanPeaks(2:end,1:end));
        % Store detections
        for k = 1:size(chanPeaks,1)
            SpikeIndex = [SpikeIndex; chanPeaks(k,2) + (CurrentBlock-1)*BlockSize*60*round(Fs)];
            ChanId = [ChanId; Chan];
            if isempty(chanPeaks(k,2:end)) 
                SpikeFV = [SpikeFV; NaN *zeros(1,8)];
            else
                SpikeFV = [SpikeFV; chanPeaks(k,2:end)];
            end
        end
    % end ChanBlock
    end
    
    
    clear EEG; clear fEEG;
end

