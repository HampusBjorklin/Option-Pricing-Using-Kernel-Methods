function [ibox,Ni] = Point2PartitionBelonging(xc,xm,Rad,Np)

% Sort point into partitions 

dim = size(xc,2);

if dim == 2
    for i = 1:Np
        flagin = sqrt((xm(i,1)-xc(:,1)).^2 + (xm(i,2)-xc(:,2)).^2) <= Rad;
        ibox(i).ind = find(flagin);
        Ni(i)=length(ibox(i).ind);
    end
    
elseif dim == 3
    for i = 1:Np
        flagin = sqrt((xm(i,1)-xc(:,1)).^2 + (xm(i,2)-xc(:,2)).^2 + (xm(i,3)-xc(:,3)).^2) <= Rad;
        ibox(i).ind = find(flagin);
        Ni(i)=length(ibox(i).ind);
    end
    
elseif dim == 4
    for i = 1:Np
        flagin = sqrt((xm(i,1)-xc(:,1)).^2 + (xm(i,2)-xc(:,2)).^2 + (xm(i,3)-xc(:,3)).^2 + (xm(i,4)-xc(:,4)).^2) <= Rad;
        ibox(i).ind = find(flagin);
        Ni(i)=length(ibox(i).ind);
    end
end