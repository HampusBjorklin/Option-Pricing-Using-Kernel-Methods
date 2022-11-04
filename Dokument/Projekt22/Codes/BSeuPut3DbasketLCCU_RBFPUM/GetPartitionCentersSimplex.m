function [xm,Np,Rad] = GetPartitionCentersSimplex(xmin,xmax,np,delta,dim)

H = (xmax-xmin)/np; % Box side
x = xmin+H/2:H:xmax-H/2; % Partition centers per dimension

if dim == 2
    [xm,ym] = ndgrid(x); % Partition grid
    xm = [xm(:) ym(:)];
    Rad = sqrt(2)*H/2*(1+delta); % Radius of partition
    ind = find(xm(:,1)+xm(:,2) <= xmax + 0.01); % Find nodes in simplex
    xm = xm(ind,:);
    
elseif dim == 3
    [xm,ym,zm] = ndgrid(x); % Partition grid
    xm = [xm(:) ym(:) zm(:)];
    Rad = sqrt(3)*H/2*(1+delta); % Radius of partition
    ind = find(xm(:,1)+xm(:,2)+xm(:,3) <= xmax + H + 0.01); % Find nodes in simplex
    xm = xm(ind,:);
    
elseif dim == 4
    [xm,ym,zm,vm] = ndgrid(x); % Partition grid
    xm = [xm(:) ym(:) zm(:) vm(:)];
    Rad = sqrt(4)*H/2*(1+delta); % Radius of partition
    ind = find(xm(:,1)+xm(:,2)+xm(:,3)+xm(:,4) <= xmax + H + 0.01); % Find nodes in simplex
    xm = xm(ind,:);
end

Np = length(xm); % # Partitions



