function [obj, gr] = obj2patchGr_R(fullfilename)

%Input: fullfilename
%�fullfilename� is optional, the user will be asked to find file if there is no input.
%Output: obj, gr

%Description: Similar to read_wobj_v2. Reads obj and creates variable obj 
%containing vertices and faces for the interface from blender. �gr� contains
%the groups made in blender in an array of structs each containing the 
%group name and associated vertices.


if ~exist('fullfilename','var')
    [file path] = uigetfile('*.obj','select .obj file');
    fullfilename = [path file];
end


fileCells = file2cellarray( fullfilename);

[ftype fdata]= fixlines(fileCells);

objects = 0;
vertices = 0;
obj.vertices = [];
which_group = 0;
gr.name = 'No Groups';
gr.vertices = [];
gr.faces = [];
g = 1;

for c = 1:size(fileCells,2)
    
    try
        type = fileCells{1,c}{1,1};
    catch
        type = 'not recognized';
    end
    

        switch type
            
            case 'g'
                
                % Checks if the the group has ended by checking if the
                % group 'name' is (null)
                if strcmp('(null)', fileCells{1,c}(1,2)) == 0 
                   
                        % Checks if the group has already been created in
                        % gr
                        if any(strcmp(fileCells{1,c}(1,2), {gr.name})) == 1
                            which_group = fileCells{1,c}{1,2};
                        else
                            gr(g).name = fileCells{1,c}(1,2);
                            gr(g).name = char(gr(g).name);
                            gr(g).vertices = [];
                            g = g + 1;
                            which_group = fileCells{1,c}{1,2};
                        end
                else
                    which_group = 0;
                end


            case 'v'

                vertices = vertices +1;

                obj.vertices(vertices,:) = [str2double(fileCells{1, c}{1,2}) str2double(fileCells{1, c}{1,3}) str2double(fileCells{1, c}{1,4})];
                
                if which_group ~= 0 
                    
                    groupPosition = find(strcmp({gr.name},{which_group}));
                    gr(groupPosition).vertices = [gr(groupPosition).vertices, str2double(fileCells{1, c}{1,2}), str2double(fileCells{1, c}{1,3}), str2double(fileCells{1, c}{1,4})];
                    
                end



            case 'o'

                objects = objects +1;

                obj.objects{objects}.type = fileCells{1, c+1}{1,1};
                
                obj.objects{objects}.vertices = [];



            case 'p'

                point = str2num(fileCells{1,c}{1,2});

                obj.objects{objects}.vertices = [obj.objects{objects}.vertices; point];


            case 'l'

                line = [str2num(fileCells{1,c}{1,2}) str2num(fileCells{1,c}{1,3})];

                obj.objects{objects}.vertices = [obj.objects{objects}.vertices; line];



            case 'f'

                face = [str2num(fileCells{1,c}{1,2}) str2num(fileCells{1,c}{1,3}) str2num(fileCells{1,c}{1,4})];

                obj.objects{objects}.vertices = [obj.objects{objects}.vertices; face];
                
                 if which_group ~= 0 
                    groupPosition = find(strcmp({gr.name},{which_group}));
                    gr(groupPosition).faces = [gr(groupPosition).faces; face];
                    face = face';
                    gr(groupPosition).vertices = cat(1, gr(groupPosition).vertices, face);
                    gr(groupPosition).vertices = unique(gr(groupPosition).vertices);
                    
                end

        end
        

    
    
    
end



















function twords=stringsplit(tline,tchar)
% Get start and end position of all "words" separated by a char
i=find(tline(2:end-1)==tchar)+1; i_start=[1 i+1]; i_end=[i-1 length(tline)];
% Create a cell array of the words
twords=cell(1,length(i_start)); for j=1:length(i_start), twords{j}=tline(i_start(j):i_end(j)); end

function file_words=file2cellarray(filename)
% Open a DI3D OBJ textfile
fid=fopen(filename,'r');
file_text=fread(fid, inf, 'uint8=>char')';
fclose(fid);
file_lines = regexp(file_text, '\n+', 'split');
file_words = regexp(file_lines, '\s+', 'split');

function [ftype fdata]=fixlines(file_words)
ftype=cell(size(file_words));
fdata=cell(size(file_words));

iline=0; jline=0;
while(iline<length(file_words))
    iline=iline+1;
    twords=removeemptycells(file_words{iline});
    if(~isempty(twords))
        % Add next line to current line when line end with '\'
        while(strcmp(twords{end},'\')&&iline<length(file_words))
            iline=iline+1;
            twords(end)=[];
            twords=[twords removeemptycells(file_words{iline})];
        end
        % Values to double
        
        type=twords{1};
        stringdold=true;
        j=0;
        switch(type)
            case{'#','$'}
                for i=2:length(twords)
                    j=j+1; twords{j}=twords{i};                    
                end    
            otherwise    
                for i=2:length(twords)
                    str=twords{i};
                    val=str2double(str);
                    stringd=~isfinite(val);
                    if(stringd)
                        j=j+1; twords{j}=str;
                    else
                        if(stringdold)
                            j=j+1; twords{j}=val;
                        else
                            twords{j}=[twords{j} val];    
                        end
                    end
                    stringdold=stringd;
                end
        end
        twords(j+1:end)=[];
        jline=jline+1;
        ftype{jline}=type;
        if(length(twords)==1), twords=twords{1}; end
        fdata{jline}=twords;
    end
end
ftype(jline+1:end)=[];
fdata(jline+1:end)=[];

function b=removeemptycells(a)
j=0; b={};
for i=1:length(a);
    if(~isempty(a{i})),j=j+1; b{j}=a{i}; end;
end

function  objects=readmtl(filename_mtl,verbose)
if(verbose),disp(['Reading Material file : ' filename_mtl]); end
file_words=file2cellarray(filename_mtl);
% Remove empty cells, merge lines split by "\" and convert strings with values to double
[ftype fdata]= fixlines(file_words);

% Surface data
objects.type(length(ftype))=0; 
objects.data(length(ftype))=0; 
no=0;
% Loop through the Wavefront object file
for iline=1:length(ftype)
  type=ftype{iline}; data=fdata{iline};
    
    % Switch on data type line
    switch(type)
        case{'#','$'}
            % Comment
            tline='  %'; 
            if(iscell(data))
                for i=1:length(data), tline=[tline ' ' data{i}]; end
            else
                tline=[tline data];
            end
            if(verbose), disp(tline); end
        case{''}
        otherwise
            no=no+1;
            if(mod(no,10000)==1), objects(no+10001).data=0; end
            objects(no).type=type;
            objects(no).data=data;
    end
end
objects=objects(1:no);
if(verbose),disp('Finished Reading Material file'); end

