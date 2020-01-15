function [sdf kernel] = msdf(psth,ftype,w,varargin)
% [sdf kernel] = msdf(psth,ftype,w,varargin)
% make spike density function
% filters a psth generated from mpsth with a Gaussian or exponential kernel to yield an SDF
% note: bin size of psth should be 1 ms!
% 
% MANDATORY ARGUMENTS
% psth      2-column vector with time base (first column) and firing rate (second column) generated from mpsth
%           alternative: single column vector, time base lacking
% ftype     'Gauss' or 'exp' for exponential or 'exGauss' / 'exGaussSimple' for exponentially modified Gaussian distribution
%           'boxcar' for a moving average window
% w         width of the filter kernel; for most purposes, reasonable values for w are 25-250 ms
%           - if ftype is 'boxcar', w denotes the width of the sliding window in ms; if w is an even number, the
%             program will instead use w+1
%           - if ftype is 'Gauss', w denotes the standard deviation (in ms), therefore 68% of the
%             spike mass will be within +- w ms
%           - if ftype is 'exp', w denotes the time constant (in ms) of the exponential distribution, i.e.
%             the time after spike onset within 63% of the mass will have decayed
%           - if ftype is 'exGauss', w should be a 2-element vector specifying a) the sd of the
%             Gaussian distribution and b) the mean of the exponential distribution (l); the higher l, the more
%             asymmetric the resulting kernel (both values to be specified in ms)
% 
% EXAMPLES
% msdf(psth,'exp',500)      returns a 2-column matrix with time base in the first column and the SDF in the second
%                           the SDF is filtered with an exponential kernel with a time constant of 500 ms
% 
% HISTORY
% March 21, 2014   included a boxcar (moving average) filter kernel and delete exGaussSimple
% Oct 30, 2013     changed output array 'sdf' in line 48 from single to double precision
% July 26, 2012    text output can be suppressed by including 'notext'
%                  w is now a two-element vector as the mean of the Gaussian was irrelevant
%                  since the code shifts the mode of the EMG distribution to zero to avoid phase shifts anyways...
%                  also outputs a reminder that the edges of the msdf are distorted because
%                  of the convolution which pads vectors with 0s
% June 28, 2012    some polishing work on the kernel output for publication
% June 8, 2012     added exGaussian
% Feb 16, 2012     added examples
% 
% by Maik C. Stüttgen, Feb 2011

%% read input argument
notext = 0;
if nargin>3
  notext = strcmp(varargin{1},'notext');
end
%% do the convolution
sdf = psth(:,end);   % preallocate for speed

switch ftype

  case 'boxcar'
    if ~mod(w,2)
      w = w + 1;
    end
    kernel      = ones(w,1)/(w);
    sdf         = conv(sdf,kernel,'same');
    [~,maxpos]  = max(kernel);
    kernel      = [(-(w-1)/2:(w-1)/2)',kernel];
  
  case 'Gauss'
    Gauss_width = max([11 6*w+1]); % hm, should be an odd number... e.g. 11
    kernel      = normpdf(-floor(Gauss_width/2):floor(Gauss_width/2),0,w);
    sdf         = conv(sdf,kernel,'same');
    [~,maxpos]  = max(kernel);
    kernel      = [-maxpos+1:length(kernel)-maxpos;kernel]';
    
    % this modification yields the same results as the previous code:
    % dummy       = conv(sdf,kernel);
    % sdf         = dummy(floor(Gauss_width/2)+1:end-floor(Gauss_width/2)); % mode of Gaussian centered on spike -> noncausal
    
  case 'exp'
    filtbase   = 1:min([3000 5*w]);
    filtfunc   = exp(-filtbase/w);
    kernel     = filtfunc./sum(filtfunc); % set integral to 1
    dummy      = conv(sdf,kernel);
    [~,maxpos] = max(kernel);
    kernel     = [-maxpos+1:length(kernel)-maxpos;kernel]';
    sdf        = dummy(1:size(psth,1));
    
  case 'exGauss'
    if numel(w)~=2
      disp('distributions not fully specified')
      return
    end
    filtbase = (0-3*w(1)):(5*w(2));
    gauss = normpdf(filtbase,0,w(1));
    expo  = exppdf(filtbase,w(2));
    kernel = conv(gauss,expo);
    [~,maxpos] = max(kernel);
    dummy  = conv(sdf,kernel);
    sdf    = dummy(maxpos:end-length(kernel)+maxpos);
    kernel = [-maxpos+1:length(kernel)-maxpos;kernel]';
    if sum(kernel(:,2))<0.99
      disp('kernel density sums to less than 0.99, adjust w')
      return
    end
  
end
%% text output
if ~notext
  disp(['--- kernel type: ' ftype ' ---'])
  pmass  = sum(kernel(1:maxpos,2));
  d5     = kernel(find(cumsum(kernel(:,2))>0.05,1,'first'),1);
  d67    = kernel(find(cumsum(kernel(:,2))>0.63,1,'first'),1);
  d90    = kernel(find(cumsum(kernel(:,2))>0.90,1,'first'),1);
  d95    = kernel(find(cumsum(kernel(:,2))>0.95,1,'first'),1);
  disp(['note that the first ' num2str(-d5,'%6.0f') ' and the last ' num2str(d95,'%6.0f')...
  ' ticks of the sdf are somewhat distorted by the filtering'])
  if nargin && strcmp(ftype,'exGauss')
    disp(['kernel has ' num2str(pmass,2) ' of its mass before 0 and rises to 0.05 until ' num2str(d5) ' ms'])
    disp(['decays to .37/.10 until ' num2str(d67) '/' num2str(d90) ' ms, respectively'])
  end
end
%% rearrange
if size(sdf,1)<size(sdf,2)
    sdf = sdf';
end
if size(psth,2)>1
  sdf = [psth(:,1) sdf];
end
