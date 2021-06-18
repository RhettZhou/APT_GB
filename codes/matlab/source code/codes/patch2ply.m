% Rhett changed 20210526, changed for blender 2.8+, not blank before \n
function patch2ply(fv,vertColors,fileName,comment)

if ~exist('fileName','var')
    [file, path] = uiputfile('*.ply','Save *.ply file to');
    
    fileName = [path file];
    
    if ~file
        disp('no file selected');
        return
    end
end



if ~exist('vertColors','var') & ~exist('fv','var')
    object = gco;
    
    fv.vertices = get(object,'vertices');
    fv.faces = get(object,'faces');
    
    
    cmap = get(gcf,'Colormap');
    lim = caxis(gca);
    
    vals = get(object,'FaceVertexCdata');


    vals(vals > lim(2)) = lim(2);
    vals(vals < lim(1)) = lim(1);
    
    vertColors = value2VertColor(vals,cmap,lim);
    
    comment = ['comment ' num2str(min(vals),2) ' to ' num2str(max(vals)) ' at/nm2']; %inserts comment about values used as limits into ply
     
end


numVerts = length(fv.vertices(:,1));
numFaces = length(fv.faces(:,1));

% indexing is 0 based in *.ply
fv.faces = fv.faces - 1;


%% writing header
header = ['ply\nformat ascii 1.0\n'];
header = [header 'comment openVA ply exporter, (c) Peter Felfer, The University of Sydney 2013, Rhett MPIE modified 2019\n'];

value = [];

if exist('comment','var')
    header = [header comment '\n'];
end

%vertex properties

header = [header 'element vertex ' num2str(numVerts) '\n'];
header = [header 'property float x\nproperty float y\nproperty float z\n'];

if exist('vertColors','var') 
   header = [header 'property uchar red\nproperty uchar green\nproperty uchar blue\n'];
    
   if max(vertColors)<=1
        header = [header 'property float value\n'];
        value = vertColors;
        vertColors = round(vertColors* 255);
    end
    
end

header = [header 'element face ' num2str(numFaces) '\n'];
header = [header 'property list uchar uint vertex_indices\n'];
header = [header 'end_header\n'];



%% file writing
if fileName
    fid = fopen(fileName,'wt');


    if( fid == -1 )
        error('Cant open file.');
        return;
    end
    
else
    return
    
end

fprintf(fid, header);


if exist('vertColors','var')
    for v = 1:numVerts
        fprintf(fid, '%f %f %f ', fv.vertices(v,:)');
        fprintf(fid, '%u %u %u', vertColors(v,:)');
        if isempty(value)
            fprintf(fid, '\n');
        else
            fprintf(fid, ' %f\n', value(v)');
        end
    end
else
    
    fprintf(fid, '%f %f %f\n', fv.vertices');
end

fprintf(fid, '3 %u %u %u\n', fv.faces');


fclose(fid);
clear fid;
end


















