function saveMarker_FFT(data, markername, fname_out)

FFT{size(data, 2)} = [];
for ipart = 1 : size(data, 2)
    try
        FFT{ipart}.(markername) = data{ipart}.(markername);
    catch
    end
end

save(fname_out, 'FFT', '-v7.3')
