function xc = GetNodesSimplex(xmin,xmax,type,N,focus,ell,mode,dim)
% Get nodes in 2-4 dimensional simplex (a corner of a cube)

if dim == 1
    switch lower(type)
        case 'uni'
            % Uniform grid with NxNxNxN points
            xc = linspace(xmin,xmax,N);
            xc = chop(xc,7);
        case 'cluster'
            % Clustered points around focus
            xc = ClusterPts(xmin,xmax,N,focus,ell,mode);
            xc = chop(xc,7);
        otherwise
            error('Requested node type not implemented')
    end   

elseif dim == 2
    switch lower(type)
        case 'uni'
            % Uniform grid with NxNxNxN points
            x = linspace(xmin,xmax,N);
            [xx,yy] = ndgrid(x);
            X = xx(:);
            Y = yy(:);
            ind_tr = find(X+Y<=xmax);
            xc = [X(ind_tr) Y(ind_tr)];
            xc = chop(xc,7);
        case 'cluster'
            % Clustered points around focus
            x = ClusterPts(xmin,xmax,N,focus,ell,mode);
            [xx,yy] = ndgrid(x);
            X = xx(:);
            Y = yy(:);
            ind_tr = find(X+Y<=xmax+0.1);
            xc = [X(ind_tr) Y(ind_tr)];
            xc = chop(xc,7);
        otherwise
            error('Requested node type not implemented')
    end
            
elseif dim == 3
    switch lower(type)
        case 'uni'
            % Uniform grid with NxNxNxN points
            x = linspace(xmin,xmax,N);
            [xx,yy,zz] = ndgrid(x);
            X = xx(:);
            Y = yy(:);
            Z = zz(:);
            ind_tr = find(X+Y+Z<=xmax);
            xc = [X(ind_tr) Y(ind_tr) Z(ind_tr)];
            xc = chop(xc,7);
        case 'cluster'
            % Clustered points around focus
            x = ClusterPts(xmin,xmax,N,focus,ell,mode);
            [xx,yy,zz] = ndgrid(x);
            X = xx(:);
            Y = yy(:);
            Z = zz(:);
            ind_tr = find(X+Y+Z<=xmax+0.1);
            xc = [X(ind_tr) Y(ind_tr) Z(ind_tr)];
            xc = chop(xc,7);
        otherwise
            error('Requested node type not implemented')
    end
            
elseif dim == 4
    switch lower(type)
        case 'uni'
            % Uniform grid with NxNxNxN points
            x = linspace(xmin,xmax,N);
            [xx,yy,zz,vv] = ndgrid(x);
            X = xx(:);
            Y = yy(:);
            Z = zz(:);
            V = vv(:);
            ind_tr = find(X+Y+Z+V<=xmax);
            xc = [X(ind_tr) Y(ind_tr) Z(ind_tr) V(ind_tr)];
            xc = round(xc,7);
        case 'cluster'
            % Clustered points around focus
            x = ClusterPts(xmin,xmax,N,focus,ell,mode);
            [xx,yy,zz,vv] = ndgrid(x);
            X = xx(:);
            Y = yy(:);
            Z = zz(:);
            V = vv(:);
            ind_tr = find(X+Y+Z+V<=xmax+0.1);
            xc = [X(ind_tr) Y(ind_tr) Z(ind_tr) V(ind_tr)];
            xc = round(xc,7);
        otherwise
            error('Requested node type not implemented')
    end
    
end
  

    
    
    
    
    
    
    
    
    
    
    
    