function [shifted, nshift, shifted_corr] = alignchannels_Xcorr(data, maxiter)

%data must be 'avg' data, ie output from ft_timelockanalysis

shifted         = data;
nshift          = [];%zeros(1, size(data.label,1));
shifted_corr    = [];%zeros(1, size(data.label,1));
toshift         = [];%zeros(1, size(data.label,1));

if size(data.label,1) == 1
  fprintf('Cannot align only one channel!\n');
  return
end


for iter = 1 : maxiter
    avg = nanmean(shifted.avg,1);
%     figure;hold;
%     plot(data.time,data.avg,'k')
%     plot(data.time,avg,'r')

    for channame = string(data.label')%ichan = 1 : size(data.label,1)
        
        nshift.(channame)       = ft_getopt(nshift, convertStringsToChars(channame), 0);
        shifted_corr.(channame) = ft_getopt(shifted_corr, convertStringsToChars(channame), 0);
        toshift.(channame)      = ft_getopt(toshift, convertStringsToChars(channame), 0);
        
        ichan = strcmp(channame, data.label);
        [X2,lags]           = nanxcorr(avg, shifted.avg(ichan,:));
        % closest peak in xcorr
        %[~, lagindx]        = findpeaks(X2);
        %[~, shiftindx]      = min(abs(lags(lagindx))); 
        %toshift.(channame)     = lags(lagindx(shiftindx));

%         % highest peak in xcorr
%         [~, lagindx]        = findpeaks(X2,'NPeaks',1,'SortStr','descend');
%         if isempty(lagindx)
%             toshift.(channame)  = NaN;
%         else
%             toshift.(channame)     = lags(lagindx);
%         end

% maximum non_sero peak
        [xcorr_allpeaks, lagindx]        = findpeaks(X2);
        
        %ignore channel if correlation is too small
        if max(xcorr_allpeaks) < 20
            toshift.(channame) = nan;
            xcorr_peak     = nan;
        else
            if isempty(lagindx)
                toshift.(channame) = nan;
                xcorr_peak     = nan;
            else
                %remove peak at zero
                xcorr_peak = xcorr_allpeaks(lags(lagindx)~=0);
                lagindx = lagindx(lags(lagindx)~=0);
                
                if isempty(lagindx)
                    toshift.(channame) = 0; %max peak is at zero
                    xcorr_peak = xcorr_allpeaks;
                else
                    [~, sel] = max(X2(lagindx));
                    lagindx = lagindx(sel);
                    xcorr_peak = xcorr_peak(sel);
                    if X2(lagindx) > X2(lags==0)/5
                        toshift.(channame) = lags(lagindx); %select a new max peak
                        xcorr_peak     = X2(lagindx);
                    else
                        toshift.(channame) = 0; %other peak is too small, keep peak at zero
                        xcorr_peak = max(xcorr_allpeaks);
                    end
                end
            end
        end
        nshift.(channame)      = nshift.(channame) + toshift.(channame);
        shifted_corr.(channame)=  xcorr_peak;
        if isnan(toshift.(channame))
            shifted.avg(ichan,:)   = nan(size(shifted.avg(ichan,:)));
        else
            shifted.avg(ichan,:)   = shift(shifted.avg(ichan,:), toshift.(channame));
        end
    end
    fprintf('Iteration %d of max %d: max absolute shift: %d samples\n', iter, maxiter, max(abs(struct2array(toshift))));
    if max(abs(struct2array(toshift))) <= 1
        fprintf('No need to continue\n');
        return
    end
end
