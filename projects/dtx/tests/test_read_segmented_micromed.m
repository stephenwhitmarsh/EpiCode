MuseStruct1      = readMuseMarkers(config{ipatient},true);
MuseStruct2 = readMuseMarkers_discontinuousMicromed(config{ipatient},true);
test1 = concatenateMuseMarker(config{ipatient},MuseStruct1,1,'SlowWave_L');
test2 = concatenateMuseMarkers_old(MuseStruct1,1,'SlowWave_L');
test3 = concatenateMuseMarker(config{ipatient},MuseStruct2,1,'SlowWave_L');
test4 = concatenateMuseMarkers_old(MuseStruct2,1,'SlowWave_L');


testdiff.clock = test2.clock-test1.clock;
test1.synctimediff = diff(test1.synctime);
test2.synctimediff = diff(test2.synctime);
max(test1.synctimediff - test2.synctimediff)
testdiff.synctime = test1.synctime-test2.synctime;
max(testdiff.synctime)
testdiff.dir = test2.dir-test1.dir;

figure;hold
plot(test1.clock,1:size(test1.clock,2))
figure
plot(test2.clock,1:size(test1.clock,2))
figure
plot(test2.synctime,1:size(test2.synctime,2))
figure
plot(test1.synctime,1:size(test1.synctime,2))