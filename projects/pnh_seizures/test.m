ntrial = length(data_all.trial);
nchans  = length(data_all.label);
spikechan = zeros(nchans,1);
for i=1:ntrial
  for j=8:8 % nchans
    spikechan(j) = spikechan(j) + all(data_all.trial{i}(j,:)==0 | data_all.trial{i}(j,:)==1 | data_all.trial{i}(j,:)==2);
    if ~all(data_all.trial{i}(j,:)==0 | data_all.trial{i}(j,:)==1 | data_all.trial{i}(j,:)==2)
        disp(i)
    end
    
    
  end
end



spikechan = (spikechan==ntrial);