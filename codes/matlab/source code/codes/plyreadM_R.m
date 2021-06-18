%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fv,fvN,v] = plyreadM_R(GBnames,pathName,smoothKey)

if ~exist('smoothKey','var')
    smoothKey = 0;
end

Count = 0;
while ~isempty(GBnames)								% split line into tokens
   [tmp,GBnames] = strtok(GBnames);
   if ~isempty(tmp)
      Count = Count + 1;                                                    % count tokens
      GB{Count} = tmp;
   end
end

fv = IEC(pathName, GB{1});

outName = [];

if smoothKey == 1
   outName = [outName 'Smooth_'];
end

outName = [outName GB{1}];

if length(GB) > 1
   for i = 2:length(GB)
      fvTemp = IEC(pathName, GB{i});
      fv = mergeStructs_R(fv,fvTemp);
      outName = [outName '_' GB{i}];
   end
end


% Find the faces from all maps. find new fv
fvAll = IEC(pathName, '0');

[VcheckAll, VcheckCur] = ismember(fvAll.vertices,fv.vertices,'rows');
VindexAll = find(VcheckAll==1);
VindexCur = VcheckCur(VcheckCur~=0);
MA2C = containers.Map(VindexAll,VindexCur);

fvCur.faces = [];
for i = 1:length(fvAll.faces)
    if ismember(fvAll.faces(i,1),VindexAll) && ismember(fvAll.faces(i,2),VindexAll) && ismember(fvAll.faces(i,3),VindexAll)
       fvCur.faces = [fvCur.faces; [MA2C(fvAll.faces(i,1)),MA2C(fvAll.faces(i,2)),MA2C(fvAll.faces(i,3))]];
    end
end

fv.faces = fvCur.faces;

if smoothKey == 1
   fv.Cs = smooth_R(fv,fv.Cs)';
   fv.IEs = smooth_R(fv,fv.IEs)';
end

v = mean(patchnormals(fv));

fCenter = mean(fv.vertices,1);

if dot(v,fCenter) < 0
   v = -v;
end

fvN = normProj(fv,v);

f1 = figure('Name','Concentration map','Visible','Off');
trisurf(fv.faces,fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3),fv.Cs);
axis equal;
rotate3d on;
shading interp;
view(v)
c = colorbar;
c.Label.String = 'Concentration (%)';
set(gca,'FontSize', 16);
fileNameOut = [pathName '\' outName '_CNorm.fig'];
set(f1, 'CreateFcn', 'set(gcbo,''Visible'',''On'')');
saveas (f1, fileNameOut);

f2 = figure('Name','Concentration map','Visible','Off');
trisurf(fvN.faces,fvN.vertices(:,1),fvN.vertices(:,2),fvN.Cs);
axis equal;
rotate3d on;
shading interp;
view([57 26])
c = colorbar;
c.Label.String = 'Concentration (%)';
set(gca,'FontSize', 16);
%hold on
%tricontour(fvN.vertices(:,1:2),fvN.faces,fvN.Cs,50);
hold off
fileNameOut = [pathName '\' outName '_C.fig'];
set(f2, 'CreateFcn', 'set(gcbo,''Visible'',''On'')');
saveas (f2, fileNameOut);

f3 = figure('Name','Interfacial excess map','Visible','Off');
trisurf(fv.faces,fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3),fv.IEs);
axis equal;
rotate3d on;
shading interp;
view(v)
c = colorbar;
c.Label.String = 'Interfacial excess (atoms/nm^{2})';
set(gca,'FontSize', 16);
fileNameOut = [pathName '\' outName '_IENorm.fig'];
set(f3, 'CreateFcn', 'set(gcbo,''Visible'',''On'')');
saveas (f3, fileNameOut);

f4 = figure('Name','Interfacial excess map','Visible','Off');
trisurf(fvN.faces,fvN.vertices(:,1),fvN.vertices(:,2),fvN.IEs);
axis equal;
rotate3d on;
shading interp;
view([57 26])
c = colorbar;
c.Label.String = 'Interfacial excess (atoms/nm^{2})';
set(gca,'FontSize', 16);
%hold on
%tricontour(fvN.vertices(:,1:2),fvN.faces,fvN.IEs,50);
hold off
fileNameOut = [pathName '\' outName '_IE.fig'];
set(f4, 'CreateFcn', 'set(gcbo,''Visible'',''On'')');
saveas (f4, fileNameOut);

end

function fv = IEC(pathName, GB)
plyN_C = [pathName '\' checkName(GB) '_C_raw.ply'];
plyN_IE = [pathName '\' checkName(GB) '_IE_raw.ply'];

if pathName(end-6) == '\'
    plyN_A = [pathName(1:end-6) 'values\' checkName(GB) '_data.txt'];
else
    plyN_A = [pathName(1:end-4) 'values\' checkName(GB) '_data.txt'];
end
%fid = fopen(plyN_A);
%tLine = fgets(fid);
fvC = geneFV(plyN_C);
fvIE = geneFV(plyN_IE);
fv.vertices = fvC.vertices;
fv.faces = fvC.faces;
fv.Cs = fvC.map;
fv.IEs = fvIE.map;
fv.As = textread(plyN_A,'%f %*f %*f %*f %*f','delimiter', ',','headerlines',0);
end


function fv = geneFV(plyN)
[Elements,Comments] = plyread(plyN);
Count = 0;
while ~isempty(Comments{1,2})								% split line into tokens
   [tmp,Comments{1,2}] = strtok(Comments{1,2});
   if ~isempty(tmp)
      Count = Count + 1;                                                    % count tokens
      Token{Count} = tmp;
   end
end

min = str2num(Token{1});
max = str2num(Token{3});
range = max - min;
map = min*ones(length(Elements.vertex.value),1) + range*Elements.vertex.value;

fv.vertices = [Elements.vertex.x,Elements.vertex.y,Elements.vertex.z];
for i = 1:length(Elements.face.vertex_indices)
    for j = 1:3
        fv.faces(i,j) = Elements.face.vertex_indices{i}(j);
    end
end
fv.faces = fv.faces + 1;
fv.map = map;
end


function name = checkName(name)
if str2num(name) <= 9 
    name = ['0' name];
end
end

function [Elements,Comments] = plyread(Path)
%PLYREAD   Read a PLY 3D data file.
%   [DATA,COMMENTS] = PLYREAD(FILENAME) reads a version 1.0 PLY file
%   FILENAME and returns a structure DATA.  The fields in this structure
%   are defined by the PLY header; each element type is a field and each
%   element property is a subfield.  If the file contains any comments,
%   they are returned in a cell string array COMMENTS.
%
%   [TRI,PTS] = PLYREAD(FILENAME,'tri') or
%   [TRI,PTS,DATA,COMMENTS] = PLYREAD(FILENAME,'tri') converts vertex
%   and face data into triangular connectivity and vertex arrays.  The
%   mesh can then be displayed using the TRISURF command.
%
%   Note: This function is slow for large mesh files (+50K faces),
%   especially when reading data with list type properties.
%
%   Example:
%   [Tri,Pts] = PLYREAD('cow.ply','tri');
%   trisurf(Tri,Pts(:,1),Pts(:,2),Pts(:,3)); 
%   colormap(gray); axis equal;
%
%   See also: PLYWRITE

% Pascal Getreuer 2004

Str = 'tri';

[fid,Msg] = fopen(Path,'rt');	% open file in read text mode

if fid == -1, error(Msg); end

Buf = fscanf(fid,'%s',1);
if ~strcmp(Buf,'ply')
   fclose(fid);
   error('Not a PLY file.'); 
end


%%% read header %%%

Position = ftell(fid);
Format = '';
NumComments = 0;
Comments = {};				% for storing any file comments
NumElements = 0;
NumProperties = 0;
Elements = [];				% structure for holding the element data
ElementCount = [];		% number of each type of element in file
PropertyTypes = [];		% corresponding structure recording property types
ElementNames = {};		% list of element names in the order they are stored in the file
PropertyNames = [];		% structure of lists of property names

while 1
   Buf = fgetl(fid);   								% read one line from file
   BufRem = Buf;
   Token = {};
   Count = 0;
   
   while ~isempty(BufRem)								% split line into tokens
      [tmp,BufRem] = strtok(BufRem);
      
      if ~isempty(tmp)
         Count = Count + 1;							% count tokens
         Token{Count} = tmp;
      end
   end
   
   if Count 		% parse line
      switch lower(Token{1})
      case 'format'		% read data format
         if Count >= 2
            Format = lower(Token{2});
            
            if Count == 3 & ~strcmp(Token{3},'1.0')
               fclose(fid);
               error('Only PLY format version 1.0 supported.');
            end
         end
      case 'comment'		% read file comment
         NumComments = NumComments + 1;
         Comments{NumComments} = '';
         for i = 2:Count
            Comments{NumComments} = [Comments{NumComments},Token{i},' '];
         end
      case 'element'		% element name
         if Count >= 3
            if isfield(Elements,Token{2})
               fclose(fid);
               error(['Duplicate element name, ''',Token{2},'''.']);
            end
            
            NumElements = NumElements + 1;
            NumProperties = 0;
   	      Elements = setfield(Elements,Token{2},[]);
            PropertyTypes = setfield(PropertyTypes,Token{2},[]);
            ElementNames{NumElements} = Token{2};
            PropertyNames = setfield(PropertyNames,Token{2},{});
            CurElement = Token{2};
            ElementCount(NumElements) = str2double(Token{3});
            
            if isnan(ElementCount(NumElements))
               fclose(fid);
               error(['Bad element definition: ',Buf]); 
            end            
         else
            error(['Bad element definition: ',Buf]);
         end         
      case 'property'	% element property
         if ~isempty(CurElement) & Count >= 3            
            NumProperties = NumProperties + 1;
            eval(['tmp=isfield(Elements.',CurElement,',Token{Count});'],...
               'fclose(fid);error([''Error reading property: '',Buf])');
            
            if tmp
               error(['Duplicate property name, ''',CurElement,'.',Token{2},'''.']);
            end            
            
            % add property subfield to Elements
            eval(['Elements.',CurElement,'.',Token{Count},'=[];'], ...
               'fclose(fid);error([''Error reading property: '',Buf])');            
            % add property subfield to PropertyTypes and save type
            eval(['PropertyTypes.',CurElement,'.',Token{Count},'={Token{2:Count-1}};'], ...
               'fclose(fid);error([''Error reading property: '',Buf])');            
            % record property name order 
            eval(['PropertyNames.',CurElement,'{NumProperties}=Token{Count};'], ...
               'fclose(fid);error([''Error reading property: '',Buf])');
         else
            fclose(fid);
            
            if isempty(CurElement)            
               error(['Property definition without element definition: ',Buf]);
            else               
               error(['Bad property definition: ',Buf]);
            end            
         end         
      case 'end_header'	% end of header, break from while loop
         break;		
      end
   end
end

%%% set reading for specified data format %%%

if isempty(Format)
	warning('Data format unspecified, assuming ASCII.');
   Format = 'ascii';
end

switch Format
case 'ascii'
   Format = 0;
case 'binary_little_endian'
   Format = 1;
case 'binary_big_endian'
   Format = 2;
otherwise
   fclose(fid);
   error(['Data format ''',Format,''' not supported.']);
end

if ~Format   
   Buf = fscanf(fid,'%f');		% read the rest of the file as ASCII data
   BufOff = 1;
else
   % reopen the file in read binary mode
   fclose(fid);
   
   if Format == 1
      fid = fopen(Path,'r','ieee-le.l64');		% little endian
   else
      fid = fopen(Path,'r','ieee-be.l64');		% big endian
   end
   
   % find the end of the header again (using ftell on the old handle doesn't give the correct position)   
   BufSize = 8192;
   Buf = [blanks(10),char(fread(fid,BufSize,'uchar')')];
   i = [];
   tmp = -11;
   
   while isempty(i)
   	i = findstr(Buf,['end_header',13,10]);			% look for end_header + CR/LF
   	i = [i,findstr(Buf,['end_header',10])];		% look for end_header + LF
      
      if isempty(i)
         tmp = tmp + BufSize;
         Buf = [Buf(BufSize+1:BufSize+10),char(fread(fid,BufSize,'uchar')')];
      end
   end
   
   % seek to just after the line feed
   fseek(fid,i + tmp + 11 + (Buf(i + 10) == 13),-1);
end


%%% read element data %%%

% PLY and MATLAB data types (for fread)
PlyTypeNames = {'char','uchar','short','ushort','int','uint','float','double', ...
   'char8','uchar8','short16','ushort16','int32','uint32','float32','double64'};
MatlabTypeNames = {'schar','uchar','int16','uint16','int32','uint32','single','double'};
SizeOf = [1,1,2,2,4,4,4,8];	% size in bytes of each type

for i = 1:NumElements
   % get current element property information
   eval(['CurPropertyNames=PropertyNames.',ElementNames{i},';']);
   eval(['CurPropertyTypes=PropertyTypes.',ElementNames{i},';']);
   NumProperties = size(CurPropertyNames,2);
   
%   fprintf('Reading %s...\n',ElementNames{i});
      
   if ~Format	%%% read ASCII data %%%
      for j = 1:NumProperties
         Token = getfield(CurPropertyTypes,CurPropertyNames{j});
         
         if strcmpi(Token{1},'list')
            Type(j) = 1;
         else
            Type(j) = 0;
			end
      end
      
      % parse buffer
      if ~any(Type)
         % no list types
         Data = reshape(Buf(BufOff:BufOff+ElementCount(i)*NumProperties-1),NumProperties,ElementCount(i))';
         BufOff = BufOff + ElementCount(i)*NumProperties;
      else
         ListData = cell(NumProperties,1);
         
         for k = 1:NumProperties
            ListData{k} = cell(ElementCount(i),1);
         end
         
         % list type
		   for j = 1:ElementCount(i)
   	      for k = 1:NumProperties
      	      if ~Type(k)
         	      Data(j,k) = Buf(BufOff);
            	   BufOff = BufOff + 1;
	            else
   	            tmp = Buf(BufOff);
      	         ListData{k}{j} = Buf(BufOff+(1:tmp))';
         	      BufOff = BufOff + tmp + 1;
            	end
            end
         end
      end
   else		%%% read binary data %%%
      % translate PLY data type names to MATLAB data type names
      ListFlag = 0;		% = 1 if there is a list type 
      SameFlag = 1;     % = 1 if all types are the same
      
      for j = 1:NumProperties
         Token = getfield(CurPropertyTypes,CurPropertyNames{j});
         
         if ~strcmp(Token{1},'list')			% non-list type
	         tmp = rem(strmatch(Token{1},PlyTypeNames,'exact')-1,8)+1;
         
            if ~isempty(tmp)
               TypeSize(j) = SizeOf(tmp);
               Type{j} = MatlabTypeNames{tmp};
               TypeSize2(j) = 0;
               Type2{j} = '';
               
               SameFlag = SameFlag & strcmp(Type{1},Type{j});
	         else
   	         fclose(fid);
               error(['Unknown property data type, ''',Token{1},''', in ', ...
                     ElementNames{i},'.',CurPropertyNames{j},'.']);
         	end
         else											% list type
            if length(Token) == 3
               ListFlag = 1;
               SameFlag = 0;
               tmp = rem(strmatch(Token{2},PlyTypeNames,'exact')-1,8)+1;
               tmp2 = rem(strmatch(Token{3},PlyTypeNames,'exact')-1,8)+1;
         
               if ~isempty(tmp) & ~isempty(tmp2)
                  TypeSize(j) = SizeOf(tmp);
                  Type{j} = MatlabTypeNames{tmp};
                  TypeSize2(j) = SizeOf(tmp2);
                  Type2{j} = MatlabTypeNames{tmp2};
	   	      else
   	   	      fclose(fid);
               	error(['Unknown property data type, ''list ',Token{2},' ',Token{3},''', in ', ...
                        ElementNames{i},'.',CurPropertyNames{j},'.']);
               end
            else
               fclose(fid);
               error(['Invalid list syntax in ',ElementNames{i},'.',CurPropertyNames{j},'.']);
            end
         end
      end
      
      % read file
      if ~ListFlag
         if SameFlag
            % no list types, all the same type (fast)
            Data = fread(fid,[NumProperties,ElementCount(i)],Type{1})';
         else
            % no list types, mixed type
            Data = zeros(ElementCount(i),NumProperties);
            
         	for j = 1:ElementCount(i)
        			for k = 1:NumProperties
               	Data(j,k) = fread(fid,1,Type{k});
              	end
         	end
         end
      else
         ListData = cell(NumProperties,1);
         
         for k = 1:NumProperties
            ListData{k} = cell(ElementCount(i),1);
         end
         
         if NumProperties == 1
            BufSize = 512;
            SkipNum = 4;
            j = 0;
            
            % list type, one property (fast if lists are usually the same length)
            while j < ElementCount(i)
               Position = ftell(fid);
               % read in BufSize count values, assuming all counts = SkipNum
               [Buf,BufSize] = fread(fid,BufSize,Type{1},SkipNum*TypeSize2(1));
               Miss = find(Buf ~= SkipNum);					% find first count that is not SkipNum
               fseek(fid,Position + TypeSize(1),-1); 		% seek back to after first count                              
               
               if isempty(Miss)									% all counts are SkipNum
                  Buf = fread(fid,[SkipNum,BufSize],[int2str(SkipNum),'*',Type2{1}],TypeSize(1))';
                  fseek(fid,-TypeSize(1),0); 				% undo last skip
                  
                  for k = 1:BufSize
                     ListData{1}{j+k} = Buf(k,:);
                  end
                  
                  j = j + BufSize;
                  BufSize = floor(1.5*BufSize);
               else
                  if Miss(1) > 1									% some counts are SkipNum
                     Buf2 = fread(fid,[SkipNum,Miss(1)-1],[int2str(SkipNum),'*',Type2{1}],TypeSize(1))';                     
                     
                     for k = 1:Miss(1)-1
                        ListData{1}{j+k} = Buf2(k,:);
                     end
                     
                     j = j + k;
                  end
                  
                  % read in the list with the missed count
                  SkipNum = Buf(Miss(1));
                  j = j + 1;
                  ListData{1}{j} = fread(fid,[1,SkipNum],Type2{1});
                  BufSize = ceil(0.6*BufSize);
               end
            end
         else
            % list type(s), multiple properties (slow)
            Data = zeros(ElementCount(i),NumProperties);
            
            for j = 1:ElementCount(i)
         		for k = 1:NumProperties
            		if isempty(Type2{k})
               		Data(j,k) = fread(fid,1,Type{k});
            		else
               		tmp = fread(fid,1,Type{k});
               		ListData{k}{j} = fread(fid,[1,tmp],Type2{k});
		            end
      		   end
      		end
         end
      end
   end
   
   % put data into Elements structure
   for k = 1:NumProperties
   	if (~Format & ~Type(k)) | (Format & isempty(Type2{k}))
      	eval(['Elements.',ElementNames{i},'.',CurPropertyNames{k},'=Data(:,k);']);
      else
      	eval(['Elements.',ElementNames{i},'.',CurPropertyNames{k},'=ListData{k};']);
		end
   end
end

clear Data ListData;
fclose(fid);

if (nargin > 1 & strcmpi(Str,'Tri')) | nargout > 2   
   % find vertex element field
   Name = {'vertex','Vertex','point','Point','pts','Pts'};
   Names = [];
   
   for i = 1:length(Name)
      if any(strcmp(ElementNames,Name{i}))
         Names = getfield(PropertyNames,Name{i});
         Name = Name{i};         
         break;
      end
   end
   
   if any(strcmp(Names,'x')) & any(strcmp(Names,'y')) & any(strcmp(Names,'z'))
      eval(['varargout{1}=[Elements.',Name,'.x,Elements.',Name,'.y,Elements.',Name,'.z];']);
   else
      varargout{1} = zeros(1,3);
	end
           
   varargout{2} = Elements;
   varargout{3} = Comments;
   Elements = [];
   
   % find face element field
   Name = {'face','Face','poly','Poly','tri','Tri'};
   Names = [];
   
   for i = 1:length(Name)
      if any(strcmp(ElementNames,Name{i}))
         Names = getfield(PropertyNames,Name{i});
         Name = Name{i};
         break;
      end
   end
   
   if ~isempty(Names)
      % find vertex indices property subfield
	   PropertyName = {'vertex_indices','vertex_indexes','vertex_index','indices','indexes'};           
      
   	for i = 1:length(PropertyName)
      	if any(strcmp(Names,PropertyName{i}))
         	PropertyName = PropertyName{i};
	         break;
   	   end
      end
      
      if ~iscell(PropertyName)
         % convert face index lists to triangular connectivity
         eval(['FaceIndices=varargout{2}.',Name,'.',PropertyName,';']);
  			N = length(FaceIndices);
   		Elements = zeros(N*2,3);
   		Extra = 0;   

			for k = 1:N
   			Elements(k,:) = FaceIndices{k}(1:3);
   
   			for j = 4:length(FaceIndices{k})
      			Extra = Extra + 1;      
	      		Elements(N + Extra,:) = [Elements(k,[1,j-1]),FaceIndices{k}(j)];
   			end
         end
         Elements = Elements(1:N+Extra,:) + 1;
      end
   end
else
   varargout{1} = Comments;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fv = normProj(fv,v)

[theta,rho,~] = cart2sph(v(1),v(2),v(3));

theta = - theta;
rho = rho - pi/2;

[fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3)] = rotZ(fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3),theta);
[fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3)] = rotY(fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3),rho);

end

function [xN,yN,zN] = rotZ(x,y,z,theta)

xN = x*cos(theta) - y*sin(theta);
yN = x*sin(theta) + y*cos(theta);
zN = z;

end

function [xN,yN,zN] = rotY(x,y,z,rho)

xN = x*cos(rho) + z*sin(rho);
yN = y;
zN = z*cos(rho) - x*sin(rho);

end
