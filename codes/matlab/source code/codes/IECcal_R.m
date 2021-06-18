function [IE, C] = IECcal_R(zPosition,cumulative,interfaceArea,interfaceLoc,zWidth,IEfigure,figKey,detectEff,atomIon)
numAtom = length(cumulative);
IEdata(:,1) = 1:numAtom;
IEdata(:,2) = cumulative;
[IEdata(:,3),IEdata(:,4),IEdata(:,5)] = ischange(cumulative,'linear','MaxNumChanges',2);          % Find the transition point
if (interfaceLoc > 1) && (interfaceLoc < numAtom)                         % for surfaces/dislocations/clusters
  if sum(IEdata(:,3)~= 0) == 2
    Lim = find(IEdata(:,3) == 1);
    Mid = floor((Lim(1) + Lim (2))/2);
    P1(1) = IEdata(1,4); P1(2) = IEdata(1,5);
    P2(1) = IEdata(numAtom,4); P2(2) = IEdata(numAtom,5);
    P0(1) = IEdata(Mid,4); P0(2) = IEdata(Mid,5);
    sl0 = abs(P2(1)-P1(1));
    sl1 = abs(P1(1)-P0(1));
    sl2 = abs(P2(1)-P0(1));
    Atoms2 = P2(1)*IEdata(Mid,1) + P2(2);
    Atoms1 = P1(1)*IEdata(Mid,1) + P1(2);
    if Atoms2 < 0 
       Atoms2 = 0;
    end
    if Atoms1 > max(cumulative)
       Atoms1 = max(cumulative);
    end
    Atoms = Atoms2 - Atoms1;
    if round((sl1+sl2)*1000) > round(sl0*1000)
      IE = (Atoms/(detectEff)) / interfaceArea; 
      zLim = find((abs(zPosition(:) - zPosition(Mid))<= zWidth/2));
      C = 100*(IEdata(max(zLim),2)-IEdata(min(zLim),2))/(IEdata(max(zLim),1)- IEdata(min(zLim),1));
      
    else
      zLim = find((abs(zPosition(:) - zPosition(interfaceLoc))<= zWidth/2));
      C = 100*(IEdata(max(zLim),2)-IEdata(min(zLim),2))/(IEdata(max(zLim),1)- IEdata(min(zLim),1));
      zLim1 = find((zPosition(:) - zPosition(1))<= zWidth);
      C1 = 100*(IEdata(max(zLim1),2)-IEdata(min(zLim1),2))/(IEdata(max(zLim1),1)- IEdata(min(zLim1),1));
      zLim2 = find((zPosition(numAtom) - zPosition(:))<= zWidth);
      C2 = 100*(IEdata(max(zLim2),2)-IEdata(min(zLim2),2))/(IEdata(max(zLim2),1)- IEdata(min(zLim2),1));
      if C ~= 0
        if C1/C < 0.5 || C1/C > 2 || C2/C < 0.5 || C2/C > 2
            IE = NaN;
        else
            IE = 0;
        end
      else
         IE = 0;       
      end
    end
    if figKey ~= 0 && ~isnan(IE)
       x1 = 0:IEdata(Lim(2),1)/100:IEdata(Lim(2),1); y1 = P1(1)*x1+P1(2);
       x2 = IEdata(Lim(1),1):(IEdata(numAtom,1)- IEdata(Lim(1),1))/100:IEdata(numAtom,1); y2 = P2(1)*x2+P2(2);
       x1 = x1/IEdata(end,1); y1 = y1/IEdata(end,2);
       x2 = x2/IEdata(end,1); y2 = y2/IEdata(end,2);
       set(gcf,'Position', [0 0 1024 768])
       subplot(2,4,figKey);
       plot(IEdata(1:numAtom,1)/IEdata(end,1), IEdata(1:numAtom,2)/IEdata(end,2),x1,y1,'--',x2,y2,'--','LineWidth',2);
       hold on
       scatter(IEdata(Mid,1)/IEdata(end,1), IEdata(Mid,2)/IEdata(end,2),'MarkerEdgeColor','y','MarkerFaceColor','r','LineWidth',2)
       hold on
       scatter(IEdata(interfaceLoc,1)/IEdata(end,1), IEdata(interfaceLoc,2)/IEdata(end,2),'MarkerEdgeColor','g','MarkerFaceColor','b','LineWidth',2)
       xlabel(['All ' atomIon 's']); ylabel(['Solute ' atomIon 's']);
       axis([0 1 0 1])
       text(0.03,0.9,['IE = ' num2str(IE,'%6.1f') ' ' atomIon 's/nm^{2}']);
       text(0.03,0.8,['C = ' num2str(C,'%4.1f') ' ' atomIon '%']);
       text(0.03,0.7,['Area = ' num2str(interfaceArea,'%4.1f') ' nm^{2}']);
       hold on
     end
  else 
     IE = NaN;
     zLim = find((abs(zPosition(:) - zPosition(interfaceLoc))<= zWidth/2));
     C = 100*(IEdata(max(zLim),2)-IEdata(min(zLim),2))/(IEdata(max(zLim),1)- IEdata(min(zLim),1));
  end
else
  IE = NaN;
  zLim = find((abs(zPosition(:) - zPosition(interfaceLoc))<= zWidth/2));
  C = 100*(IEdata(max(zLim),2)-IEdata(min(zLim),2))/(IEdata(max(zLim),1)- IEdata(min(zLim),1));
end

%%--Check: Show bad points-----------------------------------------------
aaa = 0;
if aaa == 1
    %if C < 0 && ~isnan(IE)
    if IE <= -70
        f0 = figure('Name','Ladder plots of Bad Points');
        IEdata(:,1) = 1:numAtom;
        IEdata(:,2) = cumulative;
        plot(IEdata(:,1)/IEdata(end,1),IEdata(:,2)/IEdata(end,2));
        axis equal
        xlabel(['All ' atomIon 's']); ylabel(['Solute ' atomIon 's']);
        axis([0 1 0 1])
        string1 = ['IE = ' num2str(IE,'%6.1f') ' ' atomIon 's/nm^{2}'];
        string2 = ['C = ' num2str(C,'%4.1f') ' ' atomIon '%'];
        string3 = ['Area = ' num2str(interfaceArea,'%4.1f') ' nm^{2}'];
        text(0.03, 0.9, string1);
        text(0.03, 0.8, string2);
        text(0.03, 0.7, string3);
        hold off
    end
end
end




