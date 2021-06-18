function posOut = selectRanges_R(pos,ranges)
%selects all atoms within a list of ranges of the form:
%[mcbegin mcend]
posOut = [];

for i = 1:size(ranges,1)
    multi = ranges(i,4);
    posTemp = pos(pos(:,4) > ranges(i,1) & pos(:,4) <= ranges(i,2), :);
    posTemp(:,4) = ranges(i,3);
    for j = 1:multi
     posOut = [posOut; posTemp];
    end
end