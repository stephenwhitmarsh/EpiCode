write_brainvision_eeg, lines commented with '%katia'
- modify output format from float32 to int16 to allow reading by Muse.
- add 'µV' in the header text file, otherwise Muse cannot read it


ft_writedata, line commented with '%paul'
For neuralynx data, add of this hdr field before calling write_neuralynx_ncs : 
# ncs.hdr.AcqEntName            = LABEL;
It allows Muse to get the channel name on such written data

ft_spike_maketrials, lines commented with '%paul'
Add our own fields so they are also cut in trials (samples, amplitudes)

ft_spike_select, lines commented with '%paul'
Also, add our own fields (samples, amplitudes, template, template_maxchan, cluster_group).
Besides, I changed a bit the selection of trials. Before, the trials to remove were emptied but they still existed. Here, I re-numbered the trials to they really do not exist anymore (for exemple for trials 1, 2, 3, 4, 5, if trial 4 is removed, I renumber trial 5 to trial 4, so the re-numbered trials in spike.trial are 1 2 3 4). So I can also remove trials in spike.trialinfo and spike.trialtime, which were not the case before. I did it because leaving empty trials instead of removing it could create biais in further analysis.

ft_spike_xcorrs, line commented with '%paul'
Correct one error in one specific case due to histc output if first argument is of length 1

ft_read_event : 2 lines commented with %Paul
while reading nev files, the low-function read_neuralynx_nev outputs a string associated to each event. In order to keep this information, I added it to the event structure outputed by ft_read_event