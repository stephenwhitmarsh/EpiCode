a = find(contains(t.markerlabel,'PHASE_3__START__'))
b = find(contains(t.markerlabel,'PHASE_3__END__'))
t.starttime(a(1:60))-t.starttime(b)


i = contains(t.markerlabel,'PHASE_3')

t2 = t(i,:);




i = strcmp(t2.markerlabel,'SpikeDetect')
t2(i,:) = [];
i = strcmp(t2.markerlabel,'SpikeHaT1_1')
t2(i,:) = [];
i = strcmp(t2.markerlabel,'BAD__START__')
t2(i,:) = [];
i = strcmp(t2.markerlabel,'BAD__END__')
t2(i,:) = [];
i = strcmp(t2.markerlabel,'ADStartLoss_1')
t2(i,:) = [];
i = strcmp(t2.markerlabel,'ADStartLoss_2')
t2(i,:) = [];
i = strcmp(t2.markerlabel,'ADEndLoss_1')
t2(i,:) = [];
i = strcmp(t2.markerlabel,'ADEndLoss_2')
t2(i,:) = [];
