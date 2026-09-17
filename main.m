%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for finite element calculations used in Woudenberg et al., 
% "Transgenerational polarity axis inheritance during Ceratopteris embryogenesis" 
% Nature Comm. 2026
% By Jasper van der Gucht, Physical Chemistry and Soft Matter, Wageningen
% University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Convert a B/W .bmp image into a digitalized network structure consisting of cell walls 
% Each cell wall is represented as a series of connected points, while the
% cells are triangulated to form an unstructured 2D mesh.
% This can be used as input for finite element calculations:
% The cell walls extracted from the image (called the "vertical walls") are represented as elastic beams (or springs) 
% while the triangulated cells (corresponding to the in-plane or horizontal
% cell walls) are represented as an elastic continuum in plane-stress 
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   1.  Extract cell wall network and triangulate cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Triangulations is done using DistMesh2D:
% P.-O. Persson, G. Strang, A Simple Mesh Generator in MATLAB.
% SIAM Review, Volume 46 (2), pp. 329-345, 2004
%

addpath('DistMesh2DMeshgenerator')
%%%%%%%%%%%%%%%%%%%%%%%%%

% Input is digital, thresholded microscopy image: 
%%%%%%%%%%%%%%%

imfile ='FernImage001.bmp';
networkfile = 'FernNetwork001.mat'

%Parameters for the digitalization: (see the explanation in
%"FindCellNetwork.m"):
MaxEdgeLength=5;
ExteriorIsLargestArea=true;  
MarkSpecialCells=true;
TriangulateCells=true;  
[Network,skel]= FindCellNetwork(imfile,MaxEdgeLength,ExteriorIsLargestArea,MarkSpecialCells, TriangulateCells);
%---> the output is a struct variable Network that contains all the nodes,
%segments etc (see info in FindCellNetwork.m) , and an array "skel" that is an image of the cell skeleton

plotmesh=false;
if plotmesh
    %Plot the mesh:
    X=[Network.nodes(Network.allsegments(:,1),1),Network.nodes(Network.allsegments(:,2),1)];
    Y=[Network.nodes(Network.allsegments(:,1),2),Network.nodes(Network.allsegments(:,2),2)];
    
    figure, hold on
    line(X',Y','Color','r')
    
    %Plot the segments on the boundary in blue
    Xb=[Network.nodes(Network.allsegments(Network.boundarysegments,1),1),Network.nodes(Network.allsegments(Network.boundarysegments,2),1)];
    Yb=[Network.nodes(Network.allsegments(Network.boundarysegments,1),2),Network.nodes(Network.allsegments(Network.boundarysegments,2),2)];
    line(Xb',Yb','Color','b')
end

% Save file

save(networkfile, 'Network')
disp('Image analyzed and network saved')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   2.  Finite element calculation of stress distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% PARAMETERS

par.E=1;   % Young modulus (Pa) 
par.v=0.3; % Poisson ratio
par.type=1;     % 1=plane stress; 2=plane strain 

% Parameters for balancing the contributions of vertical and horizontal
% walls
par.L=6;                % leaf thickness in micron
par.meancellsize=30;    % average cell size in horizontal direction in micron
                        % usee for scaling
par.h=0.1;              % cell wall thickness in micron




% Location of the notch (pixels)
par.rnotch0=[978,581];  
par.stiffnessgradient=true;  % cells near notch have different stiffness
% Stiffness gradient: E= par.E+ par.dE*exp[-(r/par.decaylength)^par.exp]
par.dE=4;                    % difference in stiffness between notch and far from notch
par.decaylength=130;         % decay length xi of gradient in micron
par.exp=1;                   % Shape of gradient
par.P=1e-3;                  % turgor pressure (Pa)


par.includebending=false;    %if false: walls are simple spring; otherwies beam with bending rigidity
par.includebackground=true; % if false, do not take into account the background walls

% READ GEOMETRY (NODES, EDGES, MESH)
load(networkfile)
if par.includebackground
    p=Network.allnodes;
    t=Network.triangles;        %triangular elements
else
    p=Network.nodes;
    t=[];
end
                      
N1=Network.NNodes;          %number of nodes on perpendicular walls (first N1 in p)
NT=size(p,1);               %total number of nodes
e=Network.allsegments;      %edges (beam segments)
boundary_edges=Network.boundarysegments; %edges on boundary
eb=e(boundary_edges,:);
%average cell size
areas=[];
for kk=1:Network.NCells
    areas(kk)=Network.cells{kk}.area;
end
meanarea=mean(areas);
meancellsize=sqrt(meanarea); %in pixels
par.pixelsize=par.meancellsize/meancellsize; %(micron/pixel)

p=p*par.pixelsize;
p0=p; 

par.A=par.h*par.L; %cross-sectional area of perpendicular cell walls (micron^2);
par.I=(1/12)*par.L*par.h^3;  %moment of inertia (micron^4) , needed if bending rigidity included


% INITIALIZE AND SET ESENTIAL BOUNDARIES
% Specify degrees of freedom and fixed nodes (essential boundaries)
% We will fix x,y of the lowest node, and the x-position of the highest
% node
par.rnotch=par.rnotch0;
[~,i1]=min(p(1:N1,2));
[~,i2]=max(p(1:N1,2));
if par.includebending  %in this case, the first N1 nodes also have a rotation degree of freedom
    u=[zeros(3*N1,1);zeros(2*(NT-N1),1)];
    debc=[3*(i1-1)+1,3*(i1-1)+2,3*(i2-1)+1];
    nd=3;
else
    u=zeros(2*NT,1);
    debc=[2*(i1-1)+1,2*(i1-1)+2,2*(i2-1)+1];
    nd=2;
end
Nf=size(u,1);

% STIFFNESS MATRIX
% 1. perpendicular walls
if par.includebending
    K= LinearStiffNessMatrixBeams(p,e,par,Nf);
else
    K= LinearStiffNessMatrixSprings(p,e,par,Nf); 
end

%background walls
if par.includebackground
    K2=LinearStiffnessMatrixCont(p,t,par,N1,Nf);
    K=K+K2;
end

% BOUNDARY CONDITION: PRESSURE ON OUTER WALLS
R = BoundaryStresses(p,eb,par,Nf);

% SOLVE
df=setdiff(1:Nf, debc);
du=zeros(Nf,1);
Kf = K(df, df);
Rf = R(df);
dfVals = Kf\Rf;
du(df)=dfVals;
u=u+du;

%UPDATE:
p(1:N1,1)=p(1:N1,1)+du(1:nd:nd*N1);
p(1:N1,2)=p(1:N1,2)+du(2:nd:nd*N1);
p(N1+1:end,1)=p(N1+1:end,1)+du(nd*N1+1:2:end);
p(N1+1:end,2)=p(N1+1:end,2)+du(nd*N1+2:2:end);

if par.includebending
    theta=du(3:nd:nd*N1);
else
    theta=0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   3.  Post-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate stresses and energies in elements and edges
 [results] = PostProcess(p,p0,theta,t,e,par);

% Analyze the marked (embryo) cells and store properties in cell array
% cellprops

N=length(Network.SpecialCellTriangles);
if N>0
    cellprops=cell(N);
    
    for i=1:N
        tlist=t(Network.SpecialCellTriangles{i},:);
        P1 = p0(tlist(:,1), :);
        P2 = p0(tlist(:,2), :);
        P3 = p0(tlist(:,3), :);
        
        % Triangle centroids = average of each triangle’s 3 vertices
        centroids = (P1 + P2 + P3) / 3;
        % Compute triangle areas (2D or 3D)
        P1 = [P1, zeros(size(P1,1), 1)];
        P2 = [P2, zeros(size(P2,1), 1)];
        P3 = [P3, zeros(size(P3,1), 1)];
        v1 = P2 - P1;
        v2 = P3 - P1;
        cp = cross(v1, v2, 2);  % For 2D nodes, pad to 3D with zeros if needed
        areas = 0.5 * sqrt(sum(cp.^2, 2));
        % Weighted average of centroids
        total_area = sum(areas);
        center_of_gravity = sum(centroids .* areas, 1) / total_area;

        cellprops{i}.center=center_of_gravity;
        r=par.rnotch-center_of_gravity;
        cellprops{i}.distance_to_notch=norm(r);

         %orientation with respect to notch
        theta_n = atan2(r(2), r(1));
        cellprops{i}.angle_with_respect_to_notch_direction=theta_n;
        
        %local stress state
        [el,s,h] = FindEnclosingElement(t,p0,center_of_gravity);
        Vi=results.principaldirection(el,:);
        P=results.principalstresses(el,:);
        stress_anisotropy=2*(P(1)-P(2))/(P(1)+P(2));
        theta = atan2(Vi(2), Vi(1));
        %angle with direction to notch
        dtheta = theta_n - theta;
        angle_diff = mod(dtheta + pi, 2*pi) - pi;     % range [-π, π]
        angle_deg = abs(angle_diff * 180/pi); 
        cellprops{i}.principal_stress_direction=Vi;
        cellprops{i}.principal_stress_angle=angle_deg;
        cellprops{i}.stress_anisotropy=stress_anisotropy;
        cellprops{i}.sx=results.elementstress(el,1);
        cellprops{i}.sy=results.elementstress(el,2);
        cellprops{i}.sxy=results.elementstress(el,3);
        
    end    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   3.  Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 figure, hold on
 % Define region to be plotted
 zoom=true;
 if zoom
    xmin=800;xmax=1200;
    ymin=250; ymax=700;
    
 else
    xmin=-inf;xmax=inf;
    ymin=-inf; ymax=inf;
 end
 border=10;
 % Find nodes inside the box
 X=p0(:,1);
 Y=p0(:,2);
 inside = (X >= xmin-border) & (X <= xmax+border) & (Y >= ymin-border) & (Y <= ymax+border);
 % z- value (colorcode) indicated anisotropy in  strain:
 Z=2*(results.nodalprincipalstrains(:,1)-results.nodalprincipalstrains(:,2))./(results.nodalprincipalstrains(:,1)+results.nodalprincipalstrains(:,2));

 % Keep only triangles where all 3 vertices are inside
 keep = all(inside(t), 2);  % Logical vector of triangles to keep
 % Filter TRI
 TRI_clipped = t(keep, :);
 % Plot the clipped mesh
 trisurf(TRI_clipped, X, Y, Z);

 view(2), shading interp
 %inside = (Xi >= xmin-border) & (Xi <= xmax+border) & (Yi >= ymin-border) & (Yi <= ymax+border);
    
 % plot egdes (walls):
 Xi=[p0(e(:,1),1),p0(e(:,2),1)];
 Yi=[p0(e(:,1),2),p0(e(:,2),2)];
 plot3(Xi',Yi',20*ones(size(Xi')),'Color','k','LineWidth',1)

 %Plot directions for each cell:
 V=results.principaldirection;
 d=4*par.pixelsize;  %length of line
 
 %Find cell centroids:
 centroids=[];
 for i=1:1042 %find cell centers
     x=p(Network.cells{i}.nodelist,1);
     y=p(Network.cells{i}.nodelist,2);
     if x(1) ~= x(end) || y(1) ~= y(end)
         x(end+1) = x(1);
         y(end+1) = y(1);
     end
     % Compute the signed area
     A = 0.5 * sum(x(1:end-1).*y(2:end) - x(2:end).*y(1:end-1));

     % Compute centroid coordinates
     Cx = (1/(6*A)) * sum( (x(1:end-1) + x(2:end)) .* ...
         (x(1:end-1).*y(2:end) - x(2:end).*y(1:end-1)) );
     Cy = (1/(6*A)) * sum( (y(1:end-1) + y(2:end)) .* ...
         (x(1:end-1).*y(2:end) - x(2:end).*y(1:end-1)) );

     centroids(i,:) = [Cx, Cy];
 end

 inside2=(centroids(:,1) >= xmin-border) & (centroids(:,1) <= xmax+border) & (centroids(:,2) >= ymin-border) & (centroids(:,2) <= ymax+border);
 clim([0 2])
 colorbar

 ceni=centroids(inside2,:);
 for ii=1:size(ceni,1)

     xi=ceni(ii,1);
     yi=ceni(ii,2);
     [el,s,h] = FindEnclosingElement(t,p0,[xi,yi]);
     if el>0
         Vi=V(el,:);
         xx=[xi-d*Vi(1),xi+d*Vi(1)];
         yy=[yi-d*Vi(2),yi+d*Vi(2)];
         plot3(xx,yy,20*ones(size(xx)),'k','LineWidth',1)
     end

 end
 colormap(redblue)
 if zoom
    xlim([xmin xmax])
    ylim([ymin ymax])
 end
 clim([0 0.4])




