function [fwhm1,fwhm2] = mwave(meanwave,si,varargin)
% [fwhm1,fwhm2] = mwave(meanwave,si,varargin)
% computes the widths of the first and second response peaks of
% the spike waveform (full width at half maximum)
%
% ARGUMENTS
% meanwave    required - vector containing the mean spike waveform
% si          required - scalar value indicating sampling interval in �s
%
% plot        optional - a figure will be created with the mean wave and
% pol         optional - specify polarity of the spike (direction of first signifcant deflection)
%             if not specified, mwave will choose the direction with the larger deflection
%             0: negative, 1: positive
%
% EXAMPLE
% mwave(meanwave,50,'pol',0,'plot')
%
% HISTORY
% 2014 April  added new argument 'pol'
% 2014 Jan    replaced function interp with function interp1
%
% by Maik C. St�ttgen, Feb 2013
%% upsample waveform
if ~isa(meanwave,'double');
  meanwave = double(meanwave);
end
upsamp   = 10;
dummy    = interp1(meanwave,linspace(1,length(meanwave),length(meanwave)*upsamp));
%% compute widths
% extract FWHM from upsampled waveform
p2p = max(dummy)-min(dummy); % peak-to-peak amplitude

% get polarity and plotting info
pol    = max(dummy)>abs(min(dummy));
plotit = 0;
if nargin
  for i = 1:numel(varargin)
    switch varargin{i}
      case 'pol'
        pol = varargin{i+1};
      case 'plot'
        plotit = 1;
    end
  end
end

if pol==0 % polarity negative
  % first width
  [peak1 indpeak1] = min(dummy);
  first1    = find(dummy<peak1/2,1,'first');
  last1     = find(dummy(first1+1:end)>peak1/2,1,'first')+first1;
  fwhm1     = si*(last1 - first1)/upsamp;
  % second width
  peak2     = max(dummy);
  first2    = find(dummy(indpeak1:end)>peak2/2,1,'first')+indpeak1;
  last2     = find(dummy(first2+1:end)<peak2/2,1,'first')+first2;
  try
    fwhm2   = si*(last2 - first2)/upsamp;
  end
elseif pol==1  % polarity positive
  % first width
  peak1     = max(dummy);
  first1    = find(dummy>peak1/2,1,'first');
  last1     = find(dummy(first1+1:end)<peak1/2,1,'first')+first1;
  fwhm1     = si*(last1 - first1)/upsamp;
  % second width
  peak2     = min(dummy);
  first2    = find(dummy<peak2/2,1,'first');
  last2     = find(dummy(first2+1:end)>peak2/2,1,'first')+first2;
  fwhm2     = si*(last2 - first2)/upsamp;
end
%% plot if requested
if plotit
%   figure('units','normalized','position',[.375 .375 .25 .25])
  plot(dummy),hold on
  plot([first1 first1],[-10^3 10^3],'g:')
  plot([last1 last1],[-10^3 10^3],'g:')
  plot([0 numel(dummy)],[peak1/2 peak1/2],'g:')
  plot([first1 last1],[peak1/2 peak1/2],'g')
  plot([first2 first2],[-10^3 10^3],'c:')
  plot([last2 last2],[-10^3 10^3],'c:')
  plot([0 numel(dummy)],[peak2/2 peak2/2],'c:')
  plot([first2 last2],[peak2/2 peak2/2],'c')
  axis([0 numel(dummy)+1 -(max(abs([max(dummy) min(dummy)]))) (max(abs([max(dummy) min(dummy)])))])
  set(gca,'XTick',[0,50:50:floor(numel(dummy)/50)*50],'XTickLabel',(1/upsamp)*si*[0,50:50:floor(numel(dummy)/50)*50])
  xlabel('microseconds'),ylabel('ADC units')
end
