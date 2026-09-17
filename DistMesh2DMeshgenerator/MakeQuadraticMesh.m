function [p,t] = MakeQuadraticMesh(p0,t0,nn)

% Read in triangular mesh and upgrade every triangle to a six-or 7 node
% triangle
% Places extra node on every triangle side
%   and if nn=7 also one in the middle
% p,p0: new and original list of node coordinates
% t, t0: new and original list of triangles

    Nold=size(p0,1);
    t=[t0,0*t0];  % add three columns with zeros
   
   
    t_edges= [t(:,[1,2]);
            t(:,[2,3]);
            t(:,[3,1])];
    t_edges=sort(t_edges,2);        
    [edges,ix,jx]=unique(t_edges,'rows');   %list of all unique edges in the lattice
    
    % Add nodes at midpoint of each edge:
    pnew=0.5*(p0(edges(:,1),:)+p0(edges(:,2),:));
    
    
    % Link new nodes to triangles:
    NT=size(t0,1);
    t(:,4)=Nold+jx(1:NT);
    t(:,5)=Nold+jx(NT+1:2*NT);
    t(:,6)=Nold+jx(2*NT+1:3*NT);
    p=[p0;pnew];
    
    if nn==7 
        Nn=size(p,1);
        pmid=(p0(t(:,1),:)+p0(t(:,2),:)+p0(t(:,3),:))/3;
        t(:,7)=Nn+[1:NT]';
        p=[p;pmid];
    end
