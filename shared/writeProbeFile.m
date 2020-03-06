function writeProbeFile(nb_channels, fname)

% write params file
fid = fopen(fname,'w+');

fprintf(fid,'# Probe file for Spyking Circus \n');
fprintf(fid,'# Created automatically by writeSpykingCircus.m \n \n');

fprintf(fid,'total_nb_channels = %d;\n',nb_channels);
fprintf(fid,'radius            = 10;\n\n');
fprintf(fid,'channel_groups = {\n');
fprintf(fid,'\t1: {\n');
fprintf(fid,"\t\t'channels':range(0,%d),\n",nb_channels);
fprintf(fid,"\t\t'geometry': {\n");
for i = 1 : nb_channels
    fprintf(fid,'\t\t\t\t\t\t%d: [0, %d],\n',i-1,i*50);
end
fprintf(fid,'\t\t},\n');
fprintf(fid,"\t\t'graph' : []\n");
fprintf(fid,"\t}\n");
fprintf(fid,"}\n");
fclose(fid);