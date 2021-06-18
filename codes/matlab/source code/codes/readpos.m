function [pos, varargout] = readpos(fileName)
% function [pos]=readpos(fileName)
% Reads data from a .pos file (4 floats-big endian) or .epos file
% and extract 1 matrix with x,y,z and m (Nx4)

% the filename is optional. a dialog box will pop up if no file name is
% given

if ~exist('fileName','var')
    [file path idx] = uigetfile({'*.pos';'*.epos'},'Select a pos file');
    fileName = [path file];
    disp(['file ' file ' loaded']);
end



[~ , ~, ext] = fileparts(fileName);

switch ext
    
    case '.pos'
        idx = 1;
        
    case '.POS'
        idx = 1;
        
    case '.epos'
        idx = 2;
        
    case '.EPOS'
        idx = 2;
        
    otherwise
        idx = 0;
end




if ~idx
    pos = 0;
    return
end


fid = fopen(fileName, 'r');

%% reads .pos
if idx == 1
    pos = fread(fid, inf, '4*float32', 'b');
    numAtoms = length(pos)/4;
    pos=reshape(pos, [4 numAtoms]);
    pos=pos';
    
    varargout{1} = 0;
    %% reads .epos
elseif idx == 2
    
    %% Reads through the file made of 9 floats, with 8 byte stride (the two
    %% integers at the end)
    pos = fread(fid, inf, '9*float32', 8 ,'ieee-be');
    
    
    %% reads the pulse info from epos
    fseek(fid,36,'bof');
    
    pul = fread(fid, inf, '2*uint32', 36 ,'ieee-be');

    %% Makes an array with the list of floats
    numAtoms = length(pos)/9;
    
    pos = reshape(pos, [9 numAtoms]);
    pos = single(pos);
    pos = pos';
    
    
    pulse = reshape(pul, [2,numAtoms]);
    pulse = pulse';
    
    varargout{1} = pulse;
    
    
end



fclose(fid);


end

