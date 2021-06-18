function [GBVal,interfaceAll] = IECmapGr_R(chemIn,parentPos,interfaceAll,gr,lim,rng,CVTsteps,angLim,zWidth,pathName,objName,detectEff,atomIon,allN)
% calculates an IE map for the patch 'interface' for the atoms in 'pos'
% within 'lim' nm of the interface

%%reads a pos file [x,y,z] converted to a Matlab variable and a vertex file
%%[x,y,z] and assigns every atom to the closest vertex.

%vap: m y z mc vert# disttovert distalongnormal( or to line element)
%shiftdistance
%vertex: x y z obj# nx ny nz d A(or l)

folderName = [pathName objName '\'];
if( exist(folderName, 'dir') )
   dos_cmd = sprintf( 'rmdir /S /Q "%s"',folderName);
   [~, ~] = system(dos_cmd);
end

mkdir(pathName,objName)
mkdir(pathName,[objName '\maps'])
mkdir(pathName,[objName '\values'])
mkdir(pathName,[objName '\ladder'])
mkdir(pathName,[objName '\histogram'])
mkdir(pathName,[objName '\statistics'])
mkdir(pathName,[objName '\overview'])


%% determine the ranges of solute atoms; 
% multiple ions selection is avaible; 
chem = chemIn;
Count = 0;
while ~isempty(chem)                                                       % split line into tokens, incase, we have multiple interested elements.
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

solRNG = [];
for i = 1:length(Chem)
   solRNG = [solRNG; chooseRNG_R(Chem(i),rng,atomIon)];                          % Combine all interested range file
end

%% calculate the IEmap for each group
numGr = length(gr);                                                        % number of groups;
for gri = 1:numGr                                                          % use circle to calculate the IEmap
tic
    ['Ladder plots No.' gr(gri).name]
    interface.vertices = interfaceAll.vertices(gr(gri).vertices,:);        % set interface vertices
    numGrV = length(interface.vertices);                                   % number of vertices in each group
    GrVSeq = 1:1:numGrV;                                                   % create a sequence for the vertices in each group
    mapSL = containers.Map(gr(gri).vertices,GrVSeq);                       % map vertices from each group to overall values
    mapLS = containers.Map(GrVSeq,gr(gri).vertices);                       % map vertices from overall values to each group
    numGrF = length(gr(gri).faces);                                        % number of faces in each group
    for grFi = 1:numGrF                                                    % map the faces
         grFi1 = mapSL(gr(gri).faces(grFi,1));
         grFi2 = mapSL(gr(gri).faces(grFi,2));
         grFi3 = mapSL(gr(gri).faces(grFi,3));
         interface.faces(grFi,1) = grFi1;
         interface.faces(grFi,2) = grFi2;
         interface.faces(grFi,3) = grFi3;
    end
    [areaGr,IEmapGr,CmapGr, interface] = kwikEmapMember(gr(gri).name,parentPos,interface,lim,allRNG,solRNG,CVTsteps,angLim,zWidth,pathName,objName,detectEff,atomIon);     % calculate the IEmap and Cmap in each group
    IEmapGr = fillmissing(IEmapGr,'linear','EndValues','nearest');
    CmapGr = fillmissing(CmapGr,'linear','EndValues','nearest');
    [IEAve,IEErr,CAve,CErr,AreaGB] = saveGB(gr(gri).name,areaGr,IEmapGr',CmapGr',interface,pathName,objName,atomIon);
    GBVal(:,gri) = [str2double(gr(gri).name),IEAve,IEErr,CAve,CErr,AreaGB];
    for numGrVi = 1:numGrV                                                 % map the IEmap back to the overall database
       numIEmap = mapLS(numGrVi);
       IEmap(numIEmap) = IEmapGr(numGrVi);
       Cmap(numIEmap) = CmapGr(numGrVi);
       Area(numIEmap) = areaGr(numGrVi);
       interfaceAll.vertices(numIEmap,:) = interface.vertices(numGrVi,:); % set interface vertices value back;
    end
    clear interface mapSL mapLS GrVSeq areaGr IEmapGr CmapGr
toc
end

IEmap = IEmap';
Cmap = Cmap';
Area = Area';

interfaceAll.As = Area;
interfaceAll.IEs = IEmap;
interfaceAll.Cs = Cmap;

[~,~,~,~] = saveGB('00',Area,IEmap,Cmap,interfaceAll,pathName,objName,atomIon);

%% visualising the results IEmap
f1 = figure('Name','Interfacial excess map','Visible','Off');
trisurf(interfaceAll.faces,interfaceAll.vertices(:,1),interfaceAll.vertices(:,2),interfaceAll.vertices(:,3),interfaceAll.IEs);
axis equal;
rotate3d on;
shading interp;
c = colorbar;
c.Label.String = ['Interfacial excess (' atomIon 's/nm^{2})'];
fileNameOut = [pathName objName '\overview\IEmap.fig'];
set(f1, 'CreateFcn', 'set(gcbo,''Visible'',''On'')'); 
saveas (f1, fileNameOut);

%% visualising the results Cmap
f2 = figure('Name','Grain boundary concentration map','Visible','Off');
trisurf(interfaceAll.faces,interfaceAll.vertices(:,1),interfaceAll.vertices(:,2),interfaceAll.vertices(:,3),interfaceAll.Cs);
axis equal;
rotate3d on;
shading interp;
c = colorbar;
c.Label.String = ['Concentration (' atomIon '%)'];
fileNameOut = [pathName objName '\overview\Cmap.fig'];
set(f2, 'CreateFcn', 'set(gcbo,''Visible'',''On'')');
saveas (f2, fileNameOut);

%% visualising the results
f3 = figure('Name','Interfacial excess and concentration values','Visible','Off');
yyaxis left
[~,GBValIndex] = sort(GBVal(1,:));
GBVal = GBVal(:,GBValIndex);
errorbar(GBVal(1,:)-0.1,GBVal(2,:),GBVal(3,:));
xlabel('GB Number')
ylabel(['Interfacial excess (' atomIon 's/nm^{2})'])
axis auto
xlim([min(GBVal(1,:))-1,max(GBVal(1,:))+1])
yyaxis right
errorbar(GBVal(1,:)+0.1,GBVal(4,:),GBVal(5,:));
ylabel(['Solute concentration (' atomIon '%)'])
axis auto
xlim([min(GBVal(1,:))-1,max(GBVal(1,:))+1])
hold off
fileNameOut = [pathName objName '\statistics\statistics.fig'];
set(f3, 'CreateFcn', 'set(gcbo,''Visible'',''On'')');
saveas (f3, fileNameOut);
%% Export the statistic results
IEData = [pathName '\' objName '\statistics\statistic.txt'];
csvwrite(IEData,GBVal');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [area,IEmapGr,CmapGr,interface] = kwikEmapMember(gri,parentPos,interface,lim,allRNG,solRNG,CVTsteps,angLim,zWidth,pathName,objName,detectEff,atomIon)

%% determine the parameter of CVT;

interface = centroidalVoronoiTessellation(interface,CVTsteps,angLim);

normals = patchnormals(interface);
numVerts = length(interface.vertices(:,1));

%%  volCon: volume control, cut a region that will do IEmap
volConUp = interface.vertices + lim*normals;
volConDw = interface.vertices - lim*normals;
volCon.vertices = [volConUp; volConDw];
volCon.faces = convhull(volCon.vertices); % create an alpha hull
in = inhull(parentPos(:,1:3),volCon.vertices);
posCut = parentPos(in,:);

outputPosKey = 0; % if you want to output the pos file for individual vertex region
[area,IEdata,outputPos] = vertexAreaIE_R(interface,normals,numVerts,posCut,allRNG,solRNG,outputPosKey);    

area = area';

%% calculation of IEvalue and Cmap of each individal vertex
IEFigName = ['Ladder plots No.' gri];
IEfigure = figure('Name',IEFigName,'Visible','Off');
IEFrame = 8;
IEFrameGap = ceil(numVerts/IEFrame);
figKey = 1;
figKey1 = 0;
hold on
for v = 1:numVerts
  if area(v) ~= 0
    IEdata{v} = sortrows(IEdata{v},1); 
    interfaceLoc(v) = round(median(find(abs(IEdata{v}(:,1)) == min(abs(IEdata{v}(:,1)))))); 
    cumulative{v} = cumsum(IEdata{v}(:,2));
    if mod(v,IEFrameGap) == 1 || figKey1 == 1
      [IEmapGr(v), CmapGr(v)] = IECcal_R(IEdata{v}(:,1),cumulative{v},area(v),interfaceLoc(v),zWidth,IEfigure,figKey,detectEff,atomIon);
      if isnan(IEmapGr(v)) 
         figKey1 = 1;
      else
         figKey1 = 0;
         figKey = figKey + 1;
      end
    else
      [IEmapGr(v), CmapGr(v)] = IECcal_R(IEdata{v}(:,1),cumulative{v},area(v),interfaceLoc(v),zWidth,IEfigure,0,detectEff, atomIon);
    end
  else
    if mod(v,IEFrameGap) == 1 || figKey1 == 1
       figKey1 = 1;
    end
    IEmapGr(v) = 0;
    CmapGr(v) = 0;
  end
end
hold off
fileNameOut = [pathName objName '\ladder\' IEFigName '.jpg'];
saveas (gca, fileNameOut);

%% --Find the maximum and minimum values------------------------------------

minNum = ceil(numVerts*0.01);
maxNum = ceil(numVerts*0.002);
[~,Cmax] = maxk(CmapGr,maxNum);
[~,IEmin] = mink(IEmapGr,minNum);
[~,IEmax] = maxk(IEmapGr,maxNum);
CmapGr(Cmax) = NaN;
IEmapGr(Cmax) = NaN;
IEmapGr(IEmin) = NaN;
IEmapGr(IEmax) = NaN;
CmapGr(IEmin) = NaN;
CmapGr(IEmax) = NaN;

%% --Add to calculate the IE value of bad points----------------------------
numFaces = size(interface.faces,1);
BadPCount = 0; %BadPoints = [];
CBadPCount = 0; %CBadPoints = [];
for v = 1:numVerts
 if isnan(IEmapGr(v))
    BadPContainer = []; BadPCSum = 0; BadPCSumNum = 0;
    BadPCount = BadPCount + 1;
    %BadPoints(BadPCount) = v;
    for f = 1:numFaces
       if ismember(v,interface.faces(f,:))
        BadPContainer = [BadPContainer interface.faces(f,:)];
       end
    end
    BadPContainer = unique(BadPContainer);
    BadPCCount = length (BadPContainer);
    for b = 1:BadPCCount
        if ~isnan(IEmapGr(BadPContainer(b)))
           BadPCSum = BadPCSum + IEmapGr(BadPContainer(b));
           BadPCSumNum = BadPCSumNum + 1;
        end
    end
    if BadPCSumNum ~= 0
        IEmapGr(v) = BadPCSum/BadPCSumNum;
    else
        warning('Bad Region')
    end
 end
 if isnan(CmapGr(v))
    CBadPContainer = []; CBadPCSum = 0; CBadPCSumNum = 0;
    CBadPCount = CBadPCount + 1;
    %CBadPoints(CBadPCount) = v;
    for f = 1:numFaces
       if ismember(v,interface.faces(f,:))
        CBadPContainer = [CBadPContainer interface.faces(f,:)];
       end
    end
    CBadPContainer = unique(CBadPContainer);
    CBadPCCount = length (CBadPContainer);
    for b = 1:CBadPCCount
        if ~isnan(CmapGr(CBadPContainer(b)))
           CBadPCSum = CBadPCSum + CmapGr(CBadPContainer(b));
           CBadPCSumNum = CBadPCSumNum + 1;
        end
    end
    if CBadPCSumNum ~= 0
        CmapGr(v) = CBadPCSum/CBadPCSumNum;
    else
        warning('Bad Region')
    end 
 end
end
TotalPoints = numVerts
UnresolvedIEPercentage = 100*BadPCount/numVerts
% ------------------------------------------------------------------------

end

function [IEAve,IEErr,CAve,CErr,AreaGB] = saveGB(GBID,area,IEmap,Cmap,interface,pathName,objName,atomIon)
%% export to *.ply  IEmap
AreaGB = sum(area);

minIE = min(IEmap);
maxIE = max(IEmap);
IEAve = sum(IEmap.*area)/AreaGB;
IEErr = std(IEmap,area);

limits = ['comment ' num2str(minIE,3) ' to ' num2str(maxIE,3) ' ' atomIon 's/nm2'];

IEmapPly = IEmap - minIE;
IEmapPly = IEmapPly/(maxIE-minIE);

patch2ply(interface,[IEmapPly, IEmapPly, IEmapPly],[pathName '\' objName '\maps\' GBID '_IE_raw.ply'], limits);

%% export to *.ply  Cmap
minC = min(Cmap);
maxC = max(Cmap);
CAve = sum(Cmap.*area)/AreaGB;
CErr = std(Cmap,area);

limits = ['comment ' num2str(minC,3) ' to ' num2str(maxC,3) ' ' atomIon ' percentage'];

CmapPly = Cmap - minC;
CmapPly = CmapPly/(maxC-minC);

patch2ply(interface,[CmapPly, CmapPly, CmapPly],[pathName '\' objName '\maps\' GBID '_C_raw.ply'], limits);

%% export IE and C to *.txt
IEData = [pathName '\' objName '\values\' GBID '_data.txt'];
dataMatrix = [area,IEmap,Cmap];
csvwrite(IEData,dataMatrix);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
