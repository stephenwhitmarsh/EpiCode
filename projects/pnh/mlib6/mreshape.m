function out = mreshape(data)
% function out = mreshape(data)
% 
% code accepts a one-dimensional cell array ('data') with each cell holding a vector of values
% and reshapes into a two-dimensional matrix with the first column holding all values from all vectors (across rows)
% and the second column holding a vector specifying which cell the values come from
% 
% thus, the code prepares data for analysis with kruskalwallis.m or anovan.m or mes1way.m
% 
% Maik C. Stüttgen, August 2013 @ Erasmus MC Rotterdam, The Netherlands
%% inputcheck
if ~iscell(data) || numel(size(data))>2 || all(size(data)>1)
  error ('inappropriate input format')
  return
end
%% the works
% first, count elements across cell arrays
n = zeros(numel(data)+1,1);
for i = 1:numel(data)
  n(i+1) = numel(data{i});
end
out = nan(sum(n),2); % preallocate for speed
n   = cumsum(n);
for i = 1:numel(data)
  out(n(i)+1:n(i+1),1) = reshape(data{i},[numel(data{i}),1]);
  out(n(i)+1:n(i+1),2) = i;
end