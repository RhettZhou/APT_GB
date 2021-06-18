%% Before Using Blender. by Rhett, 06/27/19
%  e-mail: x.zhou@mpie.de; modified after Peter Felfer.
%% This is a function to generate files for later imaging of atom probe data in Blender®.
%% Input parameter: 
% First variable is a string variable, i.e., ‘Au’. This is the interested ion(s), which normally decorates grain boundaries. It could be one element, ‘Au’, or multiple ions, ‘Au C2’.
% Second variable is an integer variable, i.e., 8. This indicates the percentage of atoms for imaging in Blender®.  
% Then the program requires to select one atom position file *.pos.
%  The program requires to select one mass spectrum file *.xlsx file.
%% Output parameter (stored in a subfolder “blenderImaging”); 
% A new atom position file *.pos file was generated with several percentages (i.e. 8%) randomly selected atoms restored in the new atom position file.
%  A density cloud file *Density*.raw was generated for later imaging in Blender®;
% A solute density cloud file *Solute*.raw was generated for later imaging in Blender®.
%% |Try|: S1_runMeBeforeBlender_R('Au',8);
%% 
function S1_runMeBeforeBlender_R(chemIn,percentage,pathPosName,posName,pathRngName,rngName,atomIon,allN)

chem = chemIn;

posN = [pathPosName posName];
rngN = [pathRngName rngName];

mkdir(sprintf('%s%s',pathPosName,'\blenderImaging'))

Count = 0;
while ~isempty(chem)                                                        % split line into tokens, incase, we have multiple interested elements.
   [tmp,chem] = strtok(chem);
   if ~isempty(tmp)
      Count = Count + 1;                                                   % count tokens
      Chem{Count} = tmp;
   end
end

posIn = readpos(posN);                         	                           % Input the atom position file
rng = rngRead_R(rngN);                                               % read rng file

allRNG = [];
for i = 1:length(allN)
    allRNG = [allRNG; chooseRNG_R(allN(i),rng,atomIon)];                   % Combine all interested range file
end
pos = selectRanges_R(posIn,allRNG);

solRNG = [];
for i = 1:length(Chem)
   solRNG = [solRNG; chooseRNG_R(Chem(i),rng,atomIon)];                    % Combine all interested range file
end
sol = selectRanges_R(posIn,solRNG);                                            % Crop the pos file for solute atoms

[vox, gv] = point2voxel(pos,2);                                            % Create the voxel for all atoms 
[solVox, gv] = point2voxel(sol,2,gv);                                      % Create the voxel for solute atoms

rawOutDensity = [pathPosName 'blenderImaging\' atomIon '_Density_.raw'];
rawOutSol = [pathPosName 'blenderImaging\' atomIon '_Solute_' chemIn '_.raw'];

volume2RAW(vox,gv,rawOutDensity);                                          % Generate the density file
volume2RAW(solVox,gv,rawOutSol);                                           % Generate the solute density file        

solOut = selectRanges(posIn,solRNG); 

posNew = reducePos_R(solOut,percentage);                                      % Generate reduced posfile 

fileOut = [pathPosName 'blenderImaging\' posName(1:end-4) '_' atomIon '_' chemIn '_' num2str(percentage) '%.pos'];

savepos(posNew,fileOut);

fclose('all');
