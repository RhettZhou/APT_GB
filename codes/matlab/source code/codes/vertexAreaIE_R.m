function [area,IEdata,outputPos] = vertexAreaIE_R(fv,normals,numVerts,pos,allRNG,solRNG,outputPosKey)
%calculate the area each vertex occupies. If a pos file is given, a prism
%approximation with a clipping distance clipD is used (better for atom probe data with density
%fluctuations and / or a boundary)
%outputPosKey = 1; whether output pos for individual pos
%% calcualtion of vertex area based on pos file (prism approximation)

allPos = selectRanges_R(pos,allRNG);

[~,ia,~] = unique(solRNG(:,3),'rows');
SolRNG = solRNG(ia,3);

%finding closest point for each atomic position
closest = dsearchn(fv.vertices,delaunayn(fv.vertices),allPos(:,1:3));
distVec = allPos(:,1:3) - fv.vertices(closest,:);

%distance through dot product
dist = sum(normals(closest,:) .* distVec,2);

% calculating the atomic density of the dataset

numAtom = size(allPos,1);

% calculating the area of each vertices

for vert = 1:numVerts
    vAllIn = allPos(closest == vert,1:3);
    vAllInLen = size(vAllIn,1);
    if vAllInLen > 4
       datasetVol = alphavol(vAllIn);
       deltaD = max(dist(closest == vert)) - min(dist(closest == vert));
       area(vert) = datasetVol/deltaD;
    else
       area(vert) = 0;
    end
end

numVertsCount = zeros(numVerts,1);

for i = 1:numAtom
    numVertsCount(closest(i)) = numVertsCount(closest(i)) + 1;
    IEdata{closest(i)}(numVertsCount(closest(i)),1) = dist(i);
    IEdata{closest(i)}(numVertsCount(closest(i)),2) =...
      sum(allPos(i,4)*ones(size(SolRNG,1),1) == SolRNG(:));
end
  
if outputPosKey == 1
  closestOut = dsearchn(fv.vertices,delaunayn(fv.vertices),pos(:,1:3));
  numAtomOut = size(pos,1);
  numVertsCountOut = zeros(numVerts,1);
  for i = 1:numAtomOut
      numVertsCountOut(closestOut(i)) = numVertsCountOut(closestOut(i)) + 1;
      outputPos{closestOut(i)}(numVertsCountOut(closestOut(i)),:) = pos(i,:);
  end
else
  outputPos = [];
end

end







