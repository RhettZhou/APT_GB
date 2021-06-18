function [vertColors, colBar] = value2VertColor(vals,cmap, lim)
% translates the values in vals (1 x N matrix) into vertex colors using the
% colormap 'cmap'. Default is jet. Also outputs the colorbar colorBar,
% which is a clolrmapped patch that can be used in the actual visualisation
% program.

if max(vals) < 0
    vals = uminus(vals);
end


valsN = (vals - min(lim))/(max(lim) - min(lim));


cLevels = length(cmap) -1; % 8 bit color is fully sampled


idx = uint8( round( valsN * cLevels ))+1;


vertColors = cmap(idx,:);



% creating colorbar

D = 1; %width of colorbar
H = 10; %height of colorbar


z = linspace(0,H,cLevels)';



fv.vertices = [zeros(size(z)), zeros(size(z)), z; ones(size(z))*D, zeros(size(z)), z];

fv.faces = delaunay(fv.vertices(:,1),fv.vertices(:,3));

%choice = questdlg('Export color bar to *.ply?','ply export','yes','no','no');

if false
    name = ['colorbar_' num2str(min(vals),3) '_to_' num2str(max(vals),3) '.ply'];
    [file path] = uiputfile('*.ply','save colorbar to:',name);
    patch2ply(fv,[cmap; cmap;],[path file]);
end
    


colBar.fv = fv;
colBar.colorValues = [cmap; cmap];