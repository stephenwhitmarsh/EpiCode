function saveMarker_SpikeTrials(data, markername, fname_out)

SpikeTrials{size(data, 2)} = [];
for ipart = 1 : size(data, 2)
    try
        SpikeTrials{ipart}.(markername) = data{ipart}.(markername);
    catch
    end
end

save(fname_out, 'SpikeTrials', '-v7.3')
