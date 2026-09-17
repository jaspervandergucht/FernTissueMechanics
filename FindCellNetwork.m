  
function [Network,skel]= FindCellNetwork(imfile,MaxEdgeLength,ExteriorIsLargestArea,MarkSpecialCells,TriangulateCells)


% Analyze a BW bitmap of a tissue of cell walls , and find nodes, edges, cells and joints
% Input:
%  image.bmp : filename of bmp image
%  MaxEdgeLength : maximum length of edge between two nodes
%        Inf: only corners are taken as nodes
%             otherwise extra nodes are added along cell walls if they are longer than MaxEdgeLength
%
%  MarkSpecialCells: these cells will be labeled especially by clicking
% 
%  TriangulateCells: if true, also create a background triangular mesh in
%  each cell
%
%  ExteriorIsLargestArea: if true, the code assumes that the region that
%  has the highest area is the exterior of the plant (which will not be
%  meshed); if false, it the program asks to click on it
%  
%  Output:
%
%     nodes: list of node coordinates
%     NBranchPoints: number of cornernodes; these are listed first in nodes
%
%     walls:  cell array with each cell a struct containing:
%       corners: refer to the cornernodes
%       orginalltrace: coordinates of the original wall in the skeletonized
%           image
%       onboundary: is wall on the boundary?
%       connectedcells: array with cells that share this wall (in case wall
%            is on boundary, only one cell)
%       nodelist: list of nodes along the wall [n1 n2 n3 ...nn] where n1
%            and nn are corners
%
%    walls_at_cornernode: for each node, lists connected walls
%
%    cells: cell array with each cell a struct containing
%       walls: list of walls, listed in CW-direction
%       wallorientations: for each wall, +1 is in direction of wall, -1 is in
%          oppoiste direction
%       nodelist: list of nodes along boundary (CW-direction)
%       area: area of polygon made by nodelist
%
%    joints: list of joints [n1 n2 n3] on the walls
%    branchjoints [n1, n2, n3] joints at the corners (so n2 is a
%            branchpoint)
%
%    allsegments: list of all segments [n1 n2]
%    boundarysegments: which segments in allsegments are on boundary

if nargin<5
    TriangulateCells=false;
end

if nargin<4
    MarkSpecialCells=false;
end
if nargin<3
    ExteriorIsLargestArea=true;
end

if nargin<2
   MaxEdgeLength=1; 
end


SE=strel('square',3);  
SE2=strel('disk',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. READ IMAGE AND USE WATERSHED TRANSFORM TO IDENTIFY CELLS AND EDGES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%skeleton:
I=imread(imfile);
I=im2uint8(I);   %convert to 8-bit bitmap

%invert if walls are black: (i.e. more white than black pixels)
if numel(find(I(:)==0))<numel(find(I(:)>0))
    I=~I;
end
%add black border
I = padarray(I,[3 3],0,'both');

%divide the image into cells using watershed transform:
L=watershed(I);

%Find boundaries (walls):
BW=L;
BW(BW>0)=1;     
BW=imbinarize(BW);
%BW = imclearborder(BW)
BW=~BW;  %invert: walls are 1, interiors are 0
skel=bwmorph(BW,'thin',Inf);  %skel  has 1 on walls and zeros elsewhere
L=watershed(skel,4);  %watershed again, so that boundaries of L and skel match exactly

%Find label of exterior
if ExteriorIsLargestArea
    stats=table2array(regionprops('table',L,'area'));
    [~,outerlabel]=max(stats);
else
    figure,imshow(label2rgb(L,'jet','w'));
    fprintf('Click on exterior domains that should be removed\n');
    [xi,yi] = getpts();
    for i=1:length(xi)
        outerlabel=L(round(yi(i)),round(xi(i)));
    end
    close;
end
special_label=[];
if MarkSpecialCells
    figure,imshow(label2rgb(L,'jet','w'));
    fprintf('Zoom and then press enter\n');
    zoom on;
    waitfor(gcf, 'CurrentCharacter', char(13))
    zoom reset
    zoom off
    fprintf('Click on special cells\n');
    [xi,yi] = getpts();
    for i=1:length(xi)
        special_label(i)=L(round(yi(i)),round(xi(i)));
    end
    close;
else
end

special_cells=special_label;
special_cells(special_label>outerlabel)=special_cells(special_label>outerlabel)-1;

branchpoints = bwmorph(skel, 'branchpoints');  %find branchpoints; these are the corners
endpoints = bwmorph(skel, 'endpoints');  %find endpoints

[row,col] = find(branchpoints==1);  % [row,col] coordinates of branch points
nodes=[row,col];      %row-col indices of the cornernodes
NBranchPoints=size(nodes,1);

realcells=setdiff([1:max(L(:))],outerlabel); %these are the labels of the cells, (exterior excluded)
NCells=length(realcells);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2. LINK THE BRANCH POINTS TO THE CELLS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%dilate branch points, one by one (to avoid overlapping
cells_at_nodes=cell(NBranchPoints,1);  %store which cells belong to each node
nodes_at_cell=cell(NCells+1,1);   %store which nodes belong to each cell

for i=1:NBranchPoints
    BW=zeros(size(L));
    BW(row(i),col(i))=1;  %branch point i
    %
    BW=imdilate(BW,SE);   %dilate 3x3
    conn_cells=L(BW==1);  %cells that are next to branch point
    conn_cells=unique(conn_cells(conn_cells>0));
    cells_at_nodes{i}=conn_cells';
    for j=1:length(conn_cells)
        nodes_at_cell{conn_cells(j)}=[nodes_at_cell{conn_cells(j)},i];
    end
    if length(conn_cells)<3
        fprintf('error\n');  %node must be next to at least 3 cells
    end
end

%find cells with only two nodes; these nodes have two edges between them!
nn=cellfun(@numel, nodes_at_cell);
[cell2,~]=find(nn==2);
doublenodes=[];
for i=1:length(cell2)
    doublenodes=[doublenodes; sort(nodes_at_cell{cell2(i)})];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2. CREATE LIST OF EDGES BY TRACING THE BOUNDARIES OF EACH CELL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

edges=zeros(0,2);  %store as [n1,n2] sorted by node index (n1<n2, initially)
edge_nodelist={};  %for each edge, store all nodes on the edge (incl. interior nodes)   
cell_at_edge=[];  %store cells that share this edge (max 2)
cell_edges=cell(NCells,1);  %store edges per cell, as [edge#, orientation]
ne=0;  %counter of edges

for i=1:NCells
     
    lab=realcells(i);  %label of cell i (exterior is not included here)
    
    %connected nodes:
    ni=nodes_at_cell{lab};
    %find the "lowest" node (with highest row number)
    [~,id]=max(nodes(ni,1));
    startnode=ni(id);
    %make image of cell, and dilate boundary
    bcell=zeros(size(L));
    bcell(L==lab)=1;
    for ii=1:length(ni) %(sometimes branchpoints are not exactly on boundary)
        if bcell(nodes(ni(ii),1),nodes(ni(ii),2))==1 %node is already in domain, first delete it
            bcell(nodes(ni(ii),1),nodes(ni(ii),2))=0;
        end
    end
    bcell=imdilate(bcell,SE2);  %dilate to include boundary
    
    %check whether start lies on boundary
    if bcell(nodes(startnode,1),nodes(startnode,2))==0
        % if not add point to the boundary
        bcell(nodes(startnode,1),nodes(startnode,2))=1;
    end
    %trace boundary of cell:
    Bi=bwtraceboundary(bcell,nodes(startnode,:),'W');  %trace boundary in CW direction
    %remove last if it is the same as the first
    if ismember(Bi(end,:),Bi(1,:),'rows')
        Bi(end,:)=[]; %remove last element (is same as first)
    end
    %find cornernodes on the boundary
    [found,idx]=ismember(nodes(ni,:),Bi,'rows');  %idx: index of cornernode on theboundary
    %find nodes that are not exactly on boundary (put in between two nearest points on boundary)
    tobefound=find(found==0);
    for j=1:length(tobefound)
        jj=(tobefound(j));
        nj=ni(jj);
        d=(Bi(:,1)-nodes(nj,1)).^2+(Bi(:,2)-nodes(nj,2)).^2;
        [dmin,idx2]=mink(d,2);
        idx2=sort(idx2);
        if (idx2(2)-idx2(1))==1
            %insert node and shift rest up.
            Bi=[Bi(1:idx2(1),:);
                nodes(nj,:);
                Bi(idx2(1)+1:end,:)];
            %shift indices of further nodes up
            idx(idx>idx2(1))=idx(idx>idx2(1))+1;
            idx(jj)=idx2(1)+1;
        elseif idx2(1)==1 && idx2(2)==size(Bi,1)  %insert point between first and last
            Bi=[Bi;
                nodes(nj,:)];
            idx(jj)=size(Bi,1);
        else
            fprintf('error, not two consecutive nodes in cell %d \n',i);
        end
    end
    
    %sort cornernodes along boundary, in CW-order
    [idxsorted,bb]=sort(idx);  %idsxsorted: index of the cornernodes at boundary
    ni=ni(bb);  %reorder
    
    %check whether first is indeed the start of the boundary
    %move nodes before it to the end..
    if idxsorted(1)>1
        Bi=[Bi([idxsorted(1):end],:);
            Bi([1:idxsorted-1],:)];
    end

    %split boundary up in edges:
    % for each edge store first and last node

    Nseg=length(ni);  %number of edges on cell, equal to number of corners
    i1=1;             %first node
    %edges=[];
    %edge_orientations=[];
    for j=1:Nseg
        n1=ni(i1);
        if j<Nseg   %next corner 
            i2=i1+1;
        else        %last segment connects again to first corner
            i2=1;
        end
        n2=ni(i2);

        edge_j=sort([n1,n2]);

        if n1<n2
                edge_or_j=1;
            if j<Nseg
                nlist=Bi(idxsorted(i1):idxsorted(i2),:);
            else
                nlist=[Bi(idxsorted(i1):end,:); Bi(idxsorted(1),:)];
            end
        else
            edge_or_j=-1;
            if j<Nseg
                nlist=Bi(idxsorted(i2):-1:idxsorted(i1),:);
            else
                nlist=[Bi(idxsorted(1),:);Bi(end:-1:idxsorted(i1),:)];
            end
        end
        
        ne=ne+1;
        edges(ne,:)=edge_j;  
        edge_nodelist{ne}=nlist;    
        cell_at_edge(ne)=i;  
        cell_edges{i}=[cell_edges{i}; ne,edge_or_j];  
        

        i1=i1+1;
        
    end
end

% Now duplicate edges must be found  because each edge is found twice (onve
% for each cell)

% %find unique edges,:
[edges_unique,ia,ic]=unique(edges,'rows');
[counts,ind]=histc(ic,unique(ic));
morethan2=find(counts>2);  %Edges that appear more than two times: two edges between these nodes

%Edges that appear more than two times: two edges between these nodes
for i=1:length(morethan2)
    ii=morethan2(i);
    original_edges=find(ic==ii);
    secondnode=[];
    for  i=1:length(original_edges)
        secondnode(i)=edge_nodelist{original_edges(i)}(2);
    end
    [~,e1]=find(secondnode==secondnode(1));  %edge 1
    [~,e2]=find(secondnode~=secondnode(1));  %edge 1
    ia(ii)=original_edges(e1(1));
    ic(original_edges(e1))=ii;
    %new edge:
    edges_unique=[edges_unique;edges_unique(ii,:)];
    ic(original_edges(e2))=size(edges_unique,1);
    ia(size(edges_unique,1))=original_edges(e2(1));
end

edge_nodelist_unique=edge_nodelist(ia);
%find cells attached to each edge:
cells_at_edge=cell(size(edges_unique,1),1);
for i=1:length(ic)
    cells_at_edge{ic(i)}=[cells_at_edge{ic(i)},cell_at_edge(i)];
end

for i=1:NCells
    cell_edges{i}(:,1)=ic(cell_edges{i}(:,1));   
end

nn=cellfun(@numel, cells_at_edge);  %number of cells for each edge

%find boundary edges
be=find(nn==1);
%re-orient boundary edges if needed:
for i=1:length(be)
    ci=cells_at_edge{be(i)};
    ce=cell_edges{ci};
    ind=find(ce(:,1)==be(i));
    if ce(ind,2)==-1  %re-orient edge
        cell_edges{ci}(ind,2)=1;
        edges_unique(be(i),:)=[edges_unique(be(i),2),edges_unique(be(i),1)];
        edge_nodelist_unique{be(i)}=edge_nodelist_unique{be(i)}(end:-1:1,:);
    end
end


%find all edges connected to each node
walls_at_cornernode=cell(NBranchPoints,1);
for i=1:length(edges_unique)
    walls_at_cornernode{edges_unique(i,1)}=[walls_at_cornernode{edges_unique(i,1)},i];
    walls_at_cornernode{edges_unique(i,2)}=[walls_at_cornernode{edges_unique(i,2)},i];
end

% divide walls into segements and collect data for each wall:
NWalls=length(edges_unique);
walls=cell(NWalls,1);
nn=NBranchPoints;
for i=1:length(edges_unique)
    wi.corners=edges_unique(i,:);
    wi.originaltrace=edge_nodelist_unique{i}; 
    wi.onboundary=ismember(i,be);
    wi.connectedcells=cells_at_edge{i};
    %divide wall into sub-walls
    Li=size(edge_nodelist_unique{i},1);
    nseg=ceil(Li/MaxEdgeLength); %number of segments
    t=[1:nseg-1]/nseg;
    ExtraNodesOnBranch=interparc(t,wi.originaltrace(:,1),wi.originaltrace(:,2),'linear');
    nextra=size(ExtraNodesOnBranch,1);
    nodes(nn+1:nn+nextra,:)=ExtraNodesOnBranch;
    wi.nodelist=[wi.corners(1),nn+1:nn+nextra,wi.corners(2)];
    nn=nn+nextra;
    walls{i}=wi;

end
Nnodes=size(nodes,1);
boundarywalls=be;

% collect data for each cell:
cells=cell(NCells,1);
for i=1:NCells
    cj.walls=cell_edges{i}(:,1)';
    cj.wallorientations=cell_edges{i}(:,2)';
    %nodelist:
    ni=[];
    for jj=1:length(cj.walls)
        nij=walls{cj.walls(jj)}.nodelist;
        if cj.wallorientations(jj)==1
            ni=[ni,nij(1:end-1)];
        else
            ni=[ni,nij(end:-1:2)];
        end
    end
    cj.nodelist=ni;

   cj.area=polyarea(nodes(ni,2),nodes(ni,1));
   cells{i}=cj;
end


%Find all joints:
joints=zeros(0,3);
branchjoints=zeros(0,3);
%first consider all joints in middle of walls
for i=1:NWalls
    ni=walls{i}.nodelist';
    ji=[ni(1:end-2),ni(2:end-1),ni(3:end)];
    joints=[joints;ji];
end

for i=1:NBranchPoints
    for i1=1:length(walls_at_cornernode{i})

        w1=walls{walls_at_cornernode{i}(i1)};
        ind1=find(w1.corners==i);
        if ind1==1
            n1=w1.nodelist(2);
        else
            n1=w1.nodelist(end-1);
        end
        for i2=i1+1:length(walls_at_cornernode{i})
            w2=walls{walls_at_cornernode{i}(i2)};
            ind2=find(w2.corners==i);
            if ind2==1
                n2=w2.nodelist(2);
            else
                n2=w2.nodelist(end-1);
            end
            branchjoints=[branchjoints;n1,i,n2];
        end
    end
end

%Make list of all segments
allsegments=[];
boundarysegments=[];
NSegments=0;
for i=1:NWalls
    ni=walls{i}.nodelist';
    ei=[ni(1:end-1),ni(2:end)];
    
    allsegments(NSegments+1:NSegments+size(ei,1),:)=ei;
    walls{i}.segments=[NSegments+1:NSegments+size(ei,1)];
    if walls{i}.onboundary
        boundarysegments=[boundarysegments,NSegments+1:NSegments+size(ei,1)];
    end
    NSegments=NSegments+size(ei,1);
end


%Re-orient the image, (needed because matlab treats row,col in bmp and and X,Y
% in a funy way
minx=min(nodes(:,2));
maxy=max(nodes(:,1));
p(:,1)=nodes(:,2)-minx;
p(:,2)=maxy-nodes(:,1);
Network.nodes=p;

for j=1:length(walls)
    p=[];
    p(:,1)=walls{j}.originaltrace(:,2)-minx;
    p(:,2)=maxy-walls{j}.originaltrace(:,1);
    walls{j}.originaltrace=p;

end

Network.cells=cells;
specialcelltriangles=cell(length(special_cells),1);
if TriangulateCells %loop over all cells and triangulate them
    addpath(genpath('DistMesh2D'));
    MAXIT=25;
    fh=@huniform;
    amin=MaxEdgeLength*0.9;
    allnodes=Network.nodes;
    NN=size(allnodes,1);
    triangles=[];
    for i=1:length(cells)
        
        n_i=cells{i}.nodelist;
        pfix=[Network.nodes(n_i,1),Network.nodes(n_i,2)];
        fd=@(p) dpoly(p,pfix);
        xmin=min(pfix(:,1));
        xmax=max(pfix(:,1));
        ymin=min(pfix(:,2));
        ymax=max(pfix(:,2));
        bbox=[xmin,ymin;xmax,ymax];
        if length(Network.cells{i}.walls) < length(Network.cells{i}.nodelist)
            [pi,ti]=distmesh2d_noplot(fd,fh,amin,bbox,pfix,MAXIT);
        else
            ti = delaunay(pfix(:,1), pfix(:,2));
            p1 = pfix(ti(:,1),:);
            p2 = pfix(ti(:,2),:);
            p3 = pfix(ti(:,3),:);

            area2 = (p2(:,1)-p1(:,1)).*(p3(:,2)-p1(:,2)) - ...
                    (p2(:,2)-p1(:,2)).*(p3(:,1)-p1(:,1));
            flip = area2 < 0;
            ti(flip,[2 3]) = ti(flip,[3 2]);

        end
        % find nodes in pi that match with the fixed nodes on the cell
        % boundary.
        [ind,loc] = ismember(pi,pfix,'rows');
        newnodes=pi(~ind,:);
        [tval, ~, indtval] = unique(ti);

        if length(ind) ~= length(pi)
            error('Unexpected size mismatch');
        end

        if max(tval) > length(ind)
            error('Triangulation refers to node beyond pi');
        end
        unused = setdiff(1:size(pi,1), unique(ti(:)));

        if ~isempty(unused)
            fprintf('Cell %d has %d unused mesh nodes\n',i,length(unused));
        end

        tvalnew=tval(~ind);
        newind=[1:length(tvalnew)]';
        tval(ind)=n_i(loc(tval(ind)));
        tval(~ind)=newind+NN;
        tnew=tval(indtval);
        tnew = reshape(tnew, size(ti));
        allnodes=[allnodes;newnodes];
        

        %CHEck for new boundary points, and move them away from boundary,
        %and create new element
        %find boundary nodes
        %first find edges that appear once; these are on cell boundary
        edg=[tnew(:,[1,2]);
             tnew(:,[1,3]);
             tnew(:,[2,3])];
        edg=sort(edg,2);
       [foo,ix,jx]=unique(edg,'rows');
        vec=histc(jx,1:max(jx));
        qx=find(vec==1);
        e=edg(ix(qx),:);
        %nodes on cell boundary
        b_i=unique(e(:));  
        b_ex=b_i(b_i>NN); %new nodes on boundary
        for jj=1:length(b_ex)
            %find elements that contain the node:
            ind=find(tnew(:)==b_ex(jj));
            [row,col] = ind2sub(size(tnew),ind);  %row contains elements that contain the extra node
            %find all connected nodes:
            con_n=tnew(row,:);
            con_n=con_n(:);
            con_n=sort(con_n(con_n~=b_ex(jj)));
            con_u=unique(con_n);
            ct=histc(con_n,con_u);
            con1=con_u(ct==1);  %connected once
            con2=con_u(ct>1);   %connected more than once
            %move node:
            moveto=con2(1);  %move in direction of this node
            v=allnodes(moveto,:)-allnodes(b_ex(jj),:);
            allnodes(b_ex(jj),:)=allnodes(b_ex(jj),:)+v/3;
            %create new element:
            newel=[con1(:);b_ex(jj)]';
            %check orientation:
            if simpvol(allnodes,newel)<0  %change orientation
                newel=[newel(1),newel(3),newel(2)];
            end
            tnew=[tnew;newel];
        end

        NN=size(allnodes,1);
        
        if ismember(i,special_cells)
            [~,id]=find(special_cells==i);
            specialcelltriangles{id}=[size(triangles,1)+1:size(triangles,1)+size(tnew,1)];
        end
        triangles=[triangles;tnew];

    end
    Network.triangles=triangles;
    Network.allnodes=allnodes;
    Network.NAllnodes=size(allnodes,1);
end







Network.walls=walls;
Network.cells=cells;
Network.joints=joints;
Network.branchjoints=branchjoints;
Network.allsegments=allsegments;
Network.boundarysegments=boundarysegments;
Network.walls_at_cornernode=walls_at_cornernode;
Network.NNodes=size(nodes,1);
Network.NBranchPoints=NBranchPoints;
Network.NWalls=NWalls;
Network.NCells=NCells;
Network.NSegments=size(allsegments,1);
Network.NJoints=size(joints,1);
Network.NBranchJoints=size(branchjoints,1);
Network.SpecialCells=special_cells;
Network.SpecialCellTriangles=specialcelltriangles;