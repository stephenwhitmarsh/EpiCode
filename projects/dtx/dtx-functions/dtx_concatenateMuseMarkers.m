function [markerName_Total] = dtx_concatenateMuseMarkers(cfg, MuseStruct, markerName)
%Concatenate markers of all files of the same rodent.
%Get the "clock" time of each marker position over all the files of the same rodent
%Find rodent number in dtx_setparams.m


markerName_Total = [];

Length_markerName_Total = 0;

for idir = 1:length(MuseStruct{1}) %pour chaque fichier
    if isfield(MuseStruct{1}{idir}.markers, markerName)
        if isfield(MuseStruct{1}{idir}.markers.(markerName), 'clock')
            markerName_idir = MuseStruct{1}{idir}.markers.(markerName).clock; %r�cup�rer mrk avec le m�me nom dans le dossier
            
            for iMarkerPosition = 1:length(markerName_idir)
                markerName_Total{1,iMarkerPosition+Length_markerName_Total}=markerName_idir(iMarkerPosition); %Remplir array total
                markerName_Total{2,iMarkerPosition+Length_markerName_Total}=idir;
            end
            Length_markerName_Total = Length_markerName_Total + length(markerName_idir); %se d�caler dans l'array pour remplir la suite
            
        else
            warning('Marker %s exists but is empty in %s', markerName, MuseStruct{1}{idir}.directory);
        end
    else
        warning('Marker %s does not exist in %s',markerName, MuseStruct{1}{idir}.directory);
    end
    
end
fprintf('%d times found for %s in %s \n',length(markerName_Total),markerName, cfg.prefix(1:end-1)); 