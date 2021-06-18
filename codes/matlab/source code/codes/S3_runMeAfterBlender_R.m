%% For calculating the experimental GB excesses map by Rhett, 06/19/19
%  e-mail: x.zhou@mpie.de; modified after Peter Felfer.
%% This is a function to calculate interfacial excess values and concentration within the GBs.
%% Input parameter: 
% First variable is a string variable, i.e., ‘Au’. This is the interested ion(s).
% Second variable is an integer variable, i.e., 10. This is the 'lim' nm of the interface for interfacial excess measurement.  
% Third one is used for changing the width of measured volume away from the interface, for concentration measurement, Zwidth, i.e., 1.
% Then the program requires to select an atom position file as *.pos.
% The program requires to select a mass spectrum file as *.xlsx file.
% The program requires to select a grain boundary mesh file as *.obj file.
%% Output parameter; 
% ‘Cmap.fig’, ‘IEmap.fig’, ‘statistics.fig’, and ‘statistic.txt’ are calculated maps and values for all the GBs. 
% “statistics.txt” contains five columns: GB number; IE; IE error; Concentration; Concentration error for each GB.
% The ladder diagrams of random selected 8 vertices in a GB were shown in the subfolder ‘ladder’.
% If there is a prefix ‘s_’, it is the 1st neighboring smoothed values (also averaged by its own value).
% The information of each vertex was storage in the subfolder ‘values’. Each file contains five columns: Area; IE; IE smoothed; Concentration; Concentration smoothed for each vertex;
% The patch files for each GB are storage in the folder ‘maps’. 
%% |Try|: S3_runMeAfterBlender_R('Au',10,1);
%% 
function [GBVal,Fv,SGBVal,SFv] = ...
    S3_runMeAfterBlender_R(chemIn,lim,zWidth,pathPosName,posName,pathRngName,rngName,pathObjName,objName,detectEff,atomIon,allN)

CVTsteps = 1;
angLim = 10;                                                               

posN = [pathPosName posName];
objN = [pathObjName objName];
foldName = [objName(1:end-4) '_' chemIn '_' num2str(lim) '_' num2str(zWidth) '_' atomIon];

pos = readpos(posN);
rngN = [pathRngName rngName];
rng = rngRead_R(rngN);                                                     % read rng file

[obj, gr] = obj2patchGr_R(objN);	                                       % Read the interface file
fv.vertices = obj.vertices;			                                       % Read the vertices
fv.faces = obj.objects{1}.vertices;                                        % Read the faces
grN = reGroup_R(gr);                                                       % Regroup the vertices
[GBVal,Fv] = ...
    IECmapGr_R(chemIn,pos,fv,grN,lim,rng,CVTsteps,angLim,zWidth,pathPosName,foldName,detectEff,atomIon,allN);  % Actual IE and Concentration mapping process
[SGBVal,SFv] = smoothValue_R([pathPosName foldName],atomIon);                                % Smooth the value by its surrounding, average of mean surrounding and its own value.
fclose('all');
