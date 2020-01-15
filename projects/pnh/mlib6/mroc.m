function [auroc ciLoUp] = mroc(x,y,varargin)
% function [auroc ciLoUp] = mroc(x,y,varargin)
%
% computes the area under the receiver operating characteristic (ROC) curve
% optionally, provides bootstrapped confidence intervals and plots the curve
%
% note that x and y must exclusively consist of integers!
%
% results validated a) by simulation with normcdf/norminv b) by comparison
% with web-based calculators and c) with code in the MES Toolbox
%
% INPUT ARGUMENTS
% x         input vector ('signal')
% y         input vector ('noise')
%           x and y are processed such that auroc>0.5 implies that x holds larger values and vice versa
%           so, if x is a spike count distribution during stimulation and y is a spike count distribution during baseline,
%           auroc>0.5 implies excitation and auroc <0.5 implies inhibition
%
% OPTIONAL INPUT ARGUMENTS
% plotit    if set to 1, generates a figure of the receiver operation characteristic curve
% ci        2-element vector; first element should be a scalar x 0<x<1 which determines the extent of the confidence interval;
%           second element should be an integer (recommended: 1000) specifying the number of bootstraps
%
% OUTPUT ARGUMENTS
% auroc     area under the ROC curve; equals the probability that a random sample taken from distribution x
%           is greater than a random sample taken from distribution y
% ciLoUp    2-element vector holding lower and upper bootstrapped confidence interval limits
%
% Maik C. Stüttgen, July 2012
%
% HISTORY
% Aug 2013     added optional arguments for a) plotting and b) computation of bootstrapped confidence intervals
%              added input check for integer vectors
%              added comments and documentation
%% inputcheck
% do x and y really hold integers and integers only?
if sum(x==int32(x))~=numel(x) || sum(y==int32(y))~=numel(y) || ~isvector(x) || ~isvector(y)
  error('input vectors do not exclusively consist of integers')
end
plotit = [];
ci     = [];
ciLoUp = [];
if ~isempty(varargin)
  for i = 1:2:size(varargin,2)
    switch varargin{i}
      case 'plotit'
        plotit = varargin{i+1};
      case 'ci'
        ci     = varargin{i+1};
    end
  end
end
%% the works
start_val = min([min(x) min(y)])-1;    % overall smallest value -1
end_val   = max([max(x) max(y)])-1;    % overall largest value +1
c         = start_val:end_val;         % list of criterion values
a         = zeros(numel(c),2);         % initialize HR and FA array
for i = 1:size(c,2);                   % loop through all criterion values
  a(i,1) = sum(y>c(i))/length(y);      % get p(y>c)
  a(i,2) = sum(x>c(i))/length(x);      % get p(x>c)
end
a     = [a;0,0;1,0;1,1];               % close the curve
auroc = polyarea(a(:,1),a(:,2));       % integrate
%% plot if requested
if plotit
  figure('units','normalized','position',[.4 .5 .2 .3]),hold on
  title(['AUROC=' num2str(auroc,'%1.2f')])
  scatter(a(:,1),a(:,2),'.')
  plot(a(:,1),a(:,2))
  plot([0 1],[0 1],'r')
  xlabel('p(y>c)'),ylabel('p(x>c)')
  axis([0 1 0 1])
end
%% compute bootstrapped confidence interval limits if requested
if ci
  auroc4ci = nan(ci(2),1);
  for i = 1:ci(2)
    a = zeros(numel(c),2);         % initialize HR and FA array
    x1 = randsample(x,numel(x),true);
    y1 = randsample(y,numel(y),true);
    for j = 1:size(c,2);                   % loop through all criterion values
      a(j,1) = sum(y1>c(j))/length(y1);      % get p(y>c)
      a(j,2) = sum(x1>c(j))/length(x1);      % get p(x>c)
    end
    a = [a;0,0;1,0;1,1];               % close the curve
    auroc4ci(i,1) = polyarea(a(:,1),a(:,2));       % integrate
  end
  ciLoUp = prctile(auroc4ci,100*[(1-ci(1))/2 1-(1-ci(1))/2]);
end