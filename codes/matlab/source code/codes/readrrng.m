function rng = readrrng(fileName)
%reads Cameca rrng files for analysis in Matlab as a struct

%pick file if not parsed
if ~exist('fileName','var')
    [file path idx] = uigetfile({'*.rrng';'*.RRNG'},'Select a rrng file');
    fileName = [path file];
    disp(['file ' file ' parsed']);
end

%read range file to string
fileText = fileread(fileName);

%on machines with e.g. German language systems, replace decimal ',' with
%decimal point '.'
fileText = strrep(fileText,',','.');


%find individual ranges
rangeIdx = strfind(fileText,'Range'); %find the word range in the string
rangeIdx(1) = []; %first occurance is the heading of the ranges section
numRng = length(rangeIdx);

for r = 1:numRng
    if r<numRng
        rngTxt{r} = fileText(rangeIdx(r):(rangeIdx(r+1)-1));
    else
        rngTxt{r} = fileText(rangeIdx(r):end);
    end
end

for r = 1:numRng
    rngStr = rngTxt{r};
    endHeader = strfind(rngStr,'=');
    rngStr = rngStr(endHeader+1:end);
    rngStr = strsplit(rngStr);
    
    %read range general information
    rng(r).rangeName = '';
    rng(r).mcbegin = str2num(rngStr{1});
    rng(r).mcend = str2num(rngStr{2});
    rng(r).vol = str2num(rngStr{3}(5:end)); %format: 'Vol:0.01661'
    if strcmp(rngStr{end}, '')
        rngStr(end) = [];
    end
    
    %read range ions
    ions = rngStr(4:end-1);
    numIons = length(ions);
    for i = 1:numIons
        ion = strsplit(ions{i},':');
        if strcmp(ion{1},'Name') % if range is manually named
            rng(r).atoms(i).atom = ion{2};
            rng(r).atoms(i).atomicNumber = 0;
            rng(r).atoms(i).count = 1;
            rng(r).rangeName = [ion{2} ' '];
            
            
        else
            % if range corresponds to a selected ion
            rng(r).atoms(i).atom = ion{1};
            rng(r).atoms(i).atomicNumber = symbolConvertAtomicNumber(ion{1});
            rng(r).atoms(i).count = str2num(ion{2});
            
            
            multi = rng(r).atoms(i).count;
            if multi == 1
                multi = '';
            else
                multi = num2str(multi);
            end
            
            rng(r).rangeName = [rng(r).rangeName ion{1} multi ' '];
        end
        
    end
    rng(r).rangeName(end) = [];
    
    %read range color
    color = rngStr{end}(7:end); %format: 'Color:33FFFF'
    color = [hex2dec(color(1:2)) hex2dec(color(3:4)) hex2dec(color(5:6))] / 256;
    rng(r).color = color;
    
    
    
end

%sorting by mcbegin
rng = nestedSortStruct(rng,{'mcbegin'});


end