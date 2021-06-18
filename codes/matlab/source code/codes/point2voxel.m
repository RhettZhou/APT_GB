function [vox gridVec]= point2voxel(pos,bin,gridVec)
% voxelisation of point cloud, binwidth bin, or gridvector gridVec

% only coordinates
pos = pos(:,1:3);


bin = [bin bin bin];


if ~exist('gridVec','var')
    
    
    % calculating bin centers
    binvectorRaw = linspace(0,1000*bin(1),1001);
    binvectorRaw = [fliplr(uminus(binvectorRaw(2:end))) binvectorRaw];
    
    
    %mins and maxs
    mins = min(pos) - 2*bin;
    maxs = max(pos) + 2*bin;
    
    for d = 1:3
        gridVec{d} = binvectorRaw(binvectorRaw>=mins(d) & binvectorRaw<=maxs(d));
    end
    
else
    bin = [gridVec{1}(2)-gridVec{1}(1) gridVec{2}(2)-gridVec{2}(1) gridVec{3}(2)-gridVec{3}(1)];
    
end


%% calcualting 3d histogram

% calculating bin association
for d = 1:3
    edgeVec{d} = [gridVec{d}-bin(d)/2 gridVec{d}(end)+bin(d)/2];
    [useless loc(:,d)] = histc(pos(:,d),edgeVec{d},1);
    sz(d) = length(edgeVec{d})-1;
end
clear useless

% for points on the border
sz = max([sz; max(loc,[],1)]);

% count of atoms in each voxel
vox = accumarray(loc,1,sz);


end



