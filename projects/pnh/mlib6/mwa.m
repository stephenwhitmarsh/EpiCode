function [mat opt] = mwa(tspx,tevents,trefevents,winsize,winpos,varargin)
% function [mat] = mwa(tspx,tevents,winext,winshift)
%
% code to compute the area under the ROC curve for two spike count distributions at different time points
% and different time intervals
% 
% note that function requires mroc and mnspx in Matlab path
%
% for an example, confer Figure 3 and associated manuscript text in:
% Stüttgen MC, Schwarz C (2008) Psychophysical and neurometric detection performance under stimulus uncertainty.
% Nature Neuroscience 11(9): 1091-1098.
%
% MANDATORY INPUT ARGUMENTS
% tspx            vector of spike timestamps (in seconds)
% tevents         vector of event timestamps (in seconds)
% trefevents      vector of reference event timestamps (in seconds)
% winsize         vector specifying which window sizes to use (in seconds)
% winpos          vector specifying window positions from start to end [tstart:tshift:tend],
%                 where tstart is starting time, tshift is the window displacement interval,
%                 and tend is the position of the last window (all times in seconds)
%
% OPTIONAL INPUT ARGUMENT
% plotit          if set to 1 or 2 or 3, generates a figure with AUROC values color-coded (1: grayscale, 2: autumn, 3: fancy)
%
% OUTPUT
% mat             2D-matrix with rows denoting window sizes and columns window positions (i.e. time)
% opt             combination of window size and position that yields the maximum AUROC value
%
% EXAMPLE
% [mat opt] = mwa(spx.timings,spx.eventtimings(spx.eventmarkers==9),spx.eventtimings(spx.eventmarkers==9)-2,[.05:.05:.5],[-1.5:0.05:2.5],'plotit',2);
%                 command will return a 10-by-81 matrix with 810 (approximate) AUROC values, computed from
%                 -1.5 to +2.5 seconds relative to the specified events, with window sizes ranging from 50 ms to 500 ms,
%                 and generates a fancy plot on top
%
% Maik C. Stüttgen, Summer 2013 @ Erasmus MC Rotterdam, The Netherlands
%% inputcheck
if ~isvector(tspx) || ~isvector(tevents) || ~isvector(trefevents) || ~isnumeric(winsize) || ~isnumeric(winpos)
  error('wrong input')
end
plotit   = 0;
winsize  = winsize*1000;  % because mnspx works with ms - better to scale once here than multiple times below
if ~isempty(varargin)
  for i = 1:2:size(varargin,2)
    switch varargin{i}
      case 'plotit'
        plotit = varargin{i+1};
    end
  end
end
%% the works
mat = nan(numel(winsize),numel(winpos));  % preallocate for speed
% for every winsize, compute the reference distribution (compDist)
% the compute AUROC for every winpos for the given winsize
for i = 1:numel(winsize)
  compDist = mnspx(tspx,trefevents-winsize(i)/2,winsize(i)/2,winsize(i)/2); % spike count distribution for comparison - need compute only once for each winsize
  for j = 1:numel(winpos)
    mat(i,j) = mroc(mnspx(tspx,tevents+winpos(j),winsize(i)/2,winsize(i)/2),compDist);
  end
end
%% determine optimal winsize-winpos combination
[~,Ioptwinsize] = max(max(mat,[],2));   % index of optimal window size
[~,Ioptwinpos]  = max(max(mat,[],1));   % index of optimal window position
opt = [winsize(Ioptwinsize)/1000 winpos(Ioptwinpos)];
%% plot if desired
if plotit
  figure('units','normalized','position',[.4 .5 .2 .2]),hold on
  winsize = winsize/1000;   % scale winsize back to original values for plotting
  switch plotit
    case 1
      imagesc(winpos,winsize,mat)
      colormap(gray)
    case 2
      imagesc(winpos,winsize,mat)
      colormap(autumn)
    case 3
      pcolor(winpos,winsize,mat)
      colormap(hot)
      shading interp
  end
  maxval = max(max(mat));
  minval = min(min(mat));
  r      = ceil(max([abs(0.5-maxval) abs(0.5-minval)])*10)/10;
  set(gca,'YDir','normal','Clim',[0.5-r 0.5+r])
  colorbar
  axis tight
  xlabel('window position (seconds)')
  ylabel('window size (seconds)')
  plot([0 0],[0 ceil(max(winsize))],'k-')
  plot(opt(1),opt(2),'Marker','*','MarkerEdgeColor','k')
end