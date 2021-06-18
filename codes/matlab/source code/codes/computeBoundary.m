function bnd = computeBoundary(fv)
% computes boundary vertices of a mesh
% boundary vertices are all vertices that only belong to edges that only
% belong to one triangle

DEBUG = false;

%% creating edge list
edgList = [fv.faces(:,[1 2]); fv.faces(:, [1 3]); fv.faces(:, [2 3])];
edgList = sort(edgList,2);



%% counting which ones occur once
[edgListF,b,c] = unique(edgList,'rows');

% The last column lists the counts
cnt =  accumarray(c,1);

BNDedg = edgListF(cnt == 1,:);

bnd = unique(BNDedg);


if DEBUG
    figure('Name','boundary of mesh');
    scatter3(fv.vertices(bnd,1),fv.vertices(bnd,2),fv.vertices(bnd,3));
    rotate3d on; axis equal;
end