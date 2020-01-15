function nspx = mnspx(spxtimes,trigtimes,pre,post)
% nspx = mnspx(spxtimes,trigtimes,pre,post)
% function simply returns the distribution of spike counts within a specific interval relative to a trigger
% IMPORTANT: all timestamp inputs (spxtimes, trigtimes) must be seconds and will be converted to ms in the script!
%
% MANDATORY INPUTS
% spxtimes      vector with timestamps (seconds) of spike events
% trigtimes     vector with timestamps (seconds) of trigger events
% pre           time before trigger to include in psth (milliseconds; default 1000 ms)
% post          time after trigger to include in psth (milliseconds; default 1000 ms)
%
% sorry about the two different timing formats... some compatibility issue forces me to do that... feel free to edit below!
%
% EXAMPLE
% get (and plot) the distribution of spike counts relative to the start of event 9 (food reward) and the following 2000 ms
%     nspx = mnspx(spx.timings,spx.eventtimings(spx.eventmarkers==9),0,2000)
%     hist(nspx,0.5:max(nspx)+0.5)
% 
% by Maik C. Stüttgen, Summer 2013 @ Erasmus MC Rotterdam, The Netherlands
%% preps
spxtimes  = spxtimes*1000;
trigtimes = trigtimes*1000;
nspx      = nan(numel(trigtimes),1); % preallocate for speed
%% the works
% for every trigtime, get the number of spikes in the relevant time window
for i = 1:numel(trigtimes)
  nspx(i,1) = sum(spxtimes>=trigtimes(i)-pre & spxtimes<=trigtimes(i)+post);
end