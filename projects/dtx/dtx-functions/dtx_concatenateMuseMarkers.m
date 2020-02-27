function [markerName_Total] = dtx_concatenateMuseMarkers(cfg, MuseStruct_ipart, markerName)
%Concatenate markers of all files of the same rodent.
%Only markers of one part
%Get the "clock" time of each marker position over all the files of the same rodent
%Find rodent number in dtx_setparams.m


markerName_Total = [];

Length_markerName_Total = 0;

for idir = 1:length(MuseStruct_ipart) %pour chaque fichier
    if isfield(MuseStruct_ipart{idir}.markers, markerName)
        if isfield(MuseStruct_ipart{idir}.markers.(markerName), 'clock')
            markerName_idir = MuseStruct_ipart{idir}.markers.(markerName).clock; %récupérer mrk avec le même nom dans le dossier
            
            for iMarkerPosition = 1:length(markerName_idir)
                markerName_Total{1,iMarkerPosition+Length_markerName_Total}=markerName_idir(iMarkerPosition); %Remplir array total
                markerName_Total{2,iMarkerPosition+Length_markerName_Total}=idir;
            end
            Length_markerName_Total = Length_markerName_Total + length(markerName_idir); %se décaler dans l'array pour remplir la suite
            
        else
            fprintf('Marker %s exists but is empty in %s\n', markerName, MuseStruct_ipart{idir}.directory);
        end
    else
        fprintf('Marker %s does not exist in %s\n',markerName, MuseStruct_ipart{idir}.directory);
    end
    
end
fprintf('%d times found for %s in %s \n',length(markerName_Total),markerName, cfg.prefix(1:end-1)); 
