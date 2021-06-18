function varargout = savepos(pos,fileName)

% Takes pos numatoms*4 matrix and saves it as a 
% binary file made of 4 float with big endian called fileName.pos

if ~exist('fileName','var')
    [file path ~] = uiputfile('*.pos','Save pos file');
    fileName = [path file];
    
    varargout{1} = fileName;
    
end





%% opens the file
fid = fopen(fileName,'W','b');

disp('File opened...');


%% Get the file length
numatoms=length(pos(:,1));

pos = pos';

%% reshape the matrix to get one row for posfile
mat=reshape(pos,numatoms*4,1);


%% write the matrix as floats
fwrite(fid,mat,'float32');
fclose(fid);
disp('posfile saved');

end