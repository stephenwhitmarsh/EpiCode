function [shifted, nshift] = alignXcorr(input, maxiter)

% ALIGNXCORR determines shifts according to xcorr aligment with average
%
% use as
%   [shifted, nshift] = alignXcorr(input,maxiter)
%
% This file is part of EpiCode, see
% http://www.github.com/stephenwhitmarsh/EpiCode for documentation and details.
%
%   EpiCode is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   EpiCode is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with EpiCode. If not, see <http://www.gnu.org/licenses/>.

shifted = input;
nshift  = zeros(1, size(input,1));
toshift = zeros(1, size(input,1));

if size(input,1) == 1
  fprintf('Cannot align only one trial!\n');
  return
end


for iter = 1 : maxiter
    avg = nanmean(shifted);
    for itrial = 1 : size(input,1)
        [X2,lags]           = nanxcorr(avg, shifted(itrial,:));
        [~, lagindx]        = findpeaks(X2);
        [~, shiftindx]      = min(abs(lags(lagindx))); % closest peak in xcorr
        toshift(itrial)     = lags(lagindx(shiftindx));

        % highest peak in xcorr
        % [~, lagindx]        = findpeaks(X2,'NPeaks',1,'SortStr','descend');
        % toshift(itrial)     = lags(lagindx);

        nshift(itrial)      = nshift(itrial) + toshift(itrial);
        shifted(itrial,:)   = shift(shifted(itrial,:), toshift(itrial));
    end
    fprintf('Iteration %d of max %d: max absolute shift: %d samples\n', iter, maxiter, max(abs(toshift)));
    if max(abs(toshift)) <= 1
        fprintf('No need to continue\n');
        return
    end
end