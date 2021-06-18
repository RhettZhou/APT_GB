function grN = reGroup_R(gr)

numGr = length(gr);
vertColl = [];
faceColl = [];
faceCollGr = [];
faceCollIndex = [];

for gri = 1:numGr
    grN(gri).name = gr(gri).name;
    vertColl = [vertColl; gr(gri).vertices];                               % collect all verts
    faceColl = [faceColl; gr(gri).faces];                                  % collect all faces
    faceCollGr = [faceCollGr;gri*ones(length(gr(gri).faces),1)];           % add group indexes
    a = 1:1:length(gr(gri).faces);                                         
    faceCollIndex = [faceCollIndex; a'];                                   % add face indexes in each group
end

vertColl = sort(vertColl);                                                 % sort all verts
[vertColl1,~,c] = unique(vertColl);                                        % unique sequence of verts
cntV =  accumarray(c,1);                                                   % find the frequency
GrBVert = vertColl1(cntV > 1);                                             % if verts appears multiple times

numGrBVert = length(GrBVert);                                               

GrBFaceColl = [];                                                          % face at grain boundary

for GrBVerti = 1:numGrBVert
    [fInAll,~] = find(faceColl == GrBVert(GrBVerti));                      % find all face related to the interested vert
    grVColl = faceCollGr(fInAll);
    grVColl = sort(grVColl);
    [grVColl1,~,c] = unique(grVColl);                                      % unique sequence of verts
    cntGrV =  accumarray(c,1);                                             % find the frequency
    useV = grVColl1(cntGrV == max(cntGrV));                                % maxium frequency
    if length(useV) > 1
        useV = useV(1);
    end
    for i = 1:length(fInAll)                                               % Check whether the face is a boundary face
        if faceCollGr(fInAll(i)) ~= useV
            gr(faceCollGr(fInAll(i))).faces(faceCollIndex(fInAll(i)),:) = [0 0 0];
        end
    end
end
    
for gri = 1:numGr
    grN(gri).faces = gr(gri).faces(gr(gri).faces(:,1) ~= 0,:); 
end

for gri = 1:numGr
   for i = 1:length(gr(gri).vertices)
      [key,~] = find(grN(gri).faces == gr(gri).vertices(i));
      if isempty(key)
          gr(gri).vertices(i) = 0;
      end
   end
end

for gri = 1:numGr
   grN(gri).vertices = gr(gri).vertices(gr(gri).vertices ~= 0); 
end

end