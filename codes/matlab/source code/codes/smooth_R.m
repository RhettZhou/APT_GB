function fN = smooth_R(fv,f)
numV = length(fv.vertices);
for Vi = 1:numV
    Vcoll = [];
    [fVi,~] = find(fv.faces == Vi);                                        % find all face related to the interested vert
    for i = 1:3
        Vcoll = [Vcoll; fv.faces(fVi,i)];
    end
    Vcoll = sort(Vcoll);
    [Vcoll1,~,~] = unique(Vcoll);
    fN(Vi) = (mean(f(Vcoll1))+f(Vi))/2;                                    % Averaged too much, reduce the effect of smooth
end



