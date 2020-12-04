function [X2,lags] = nanxcorr(data1,varargin)


scaled1 = (data1 - nanmean(data1)) / nanstd(data1);
scaled1(isnan(data1)) = 0;

if nargin > 1
    if isscalar(varargin{1})
        [X2,lags] = xcorr(scaled1,varargin{1});
    else
        data2 = varargin{1};
        scaled2 = (data2 - nanmean(data2)) / nanstd(data2);
        scaled2(isnan(data2)) = 0;
        [X2,lags] = xcorr(scaled1,scaled2);
    end
end

if nargin > 2
    data2 = varargin{1};
    nlags = varargin{2};
    scaled2 = (data2 - nanmean(data2)) / nanstd(data2);
    scaled2(isnan(data2)) = 0;
    [X2,lags] = xcorr(scaled1,scaled2,nlags);
end



% 
% 
% if nargin > 1
%     
%     [X2,lags] = xcorr(scaled1,scaled2,nlags);
%     e
%     [X2,lags] = xcorr(scaled1,scaled2);
%     
%     
    % out=zeros(nlags,1);
    % out(1)=1;
    % for i=2:nlags+1
    %     out(i)=corr(data(i:end),data(1:end-i+1),'rows','complete');
    % end
    % stem(0:nlags,out)
    % title('sample ACF')