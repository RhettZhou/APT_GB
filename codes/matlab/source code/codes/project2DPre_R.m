%% Find the idea theta value. by Rhett, 09/17/19
%  e-mail: x.zhou@mpie.de;
%% This is a function to find find the C information of the tip
%% 
function [pAll,pSol,meanAll,stdAll] = ...
    project2DPre_R(chemIn,percentage,resolu,pathPosName,posName,pathRngName,rngName,atomIon,allN)

posN = [pathPosName posName];
rngN = [pathRngName rngName];
rng = rngRead_R(rngN);

pos = readpos(posN);                         	                           % Input the atom position file
posRedu = reducePos_R(pos,percentage);                                     % Generate reduced posfile 

chem = chemIn;
Count = 0;
while ~isempty(chem)                                                        % split line into tokens
   [tmp,chem] = strtok(chem);
   if ~isempty(tmp)
      Count = Count + 1;                                                   % count tokens
      Chem{Count} = tmp;
   end
end

allRNG = [];
for i = 1:length(allN)
    allRNG = [allRNG; chooseRNG_R(allN(i),rng,atomIon)];                   % Combine all interested range file
end
pAll = selectRanges_R(posRedu,allRNG);

solRNG = [];
for i = 1:length(Chem)
   solRNG = [solRNG; chooseRNG_R(Chem(i),rng,atomIon)];                          % Combine all interested range file
end
pSol = selectRanges_R(posRedu,solRNG);

pAll = pAll(:,1:3);
pSol = pSol(:,1:3);

[vox, gv] = point2voxel(pAll,5*resolu);                                    % Create the voxel for all atoms 
[solVox, ~] = point2voxel(pSol,5*resolu,gv);                               % Create the voxel for solute atoms
cVox = solVox./vox;
stdAll = std(cVox(:), 'omitnan');
meanAll = mean(cVox(:), 'omitnan');

end