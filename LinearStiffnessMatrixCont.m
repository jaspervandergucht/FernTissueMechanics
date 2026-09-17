function [K] = LinearStiffnessMatrixCont(p,t,par,N1,Nf)

% Calculate the linear stiffness matrix for 3node  triangular mesh
% for plane stress (stress_type=1) or plane strain (2) 
% KL=int(BL'*C*BL)dV=element stiffness matrix
% uses Gauss integration and transformation to local coordinates
%



if nargin<5
    Nf=2*size(p,1);
end
if nargin<4
    N1=Nf/2;
end
Nt=size(t,1);


switch (par.type)  %1. plane stress; 2. plane strain
case 1
    C = par.E/(1 - par.v^2)*[1, par.v, 0; par.v, 1, 0; 0, 0, (1 - par.v)/2];
case 2
    C = par.E/((1 + par.v)*(1 - 2*par.v))*[1 - par.v, par.v, 0; par.v, 1 - par.v, 0;
        0, 0, (1 - 2*par.v)/2]; 
end  

% define integration points for Gauss quadrature
% int F ds dt = 0.5*sum w_i*F(ri,si)
% 3 node triangles --> 1 GP

gpLocs = [1/3,1/3];  
gpWts = [1];

%%%%%%%%%%%%

% matrix K is stored as three columns [I,J,Val]: row, col, value of non-zero elememts
% http://blogs.mathworks.com/loren/2007/03/01/creating-sparse-finite-element-matrices-in-matlab/


KI = zeros (10, 1) ;    %set initial size to 10 (grow later when needed)
KJ = zeros (10, 1) ;    
KX = zeros (10, 1) ;
ntriplets = 0 ;        %number of entries

if par.includebending
    nd=3;
else
    nd=2;
end
% Loop over all elements
for e=1:Nt
    lm=t(e,:);  %node indices
    kl=zeros(6);
    lmg=[];
    for j=1:3
        if lm(j)<N1+1 
            lmg=[lmg,nd*(lm(j)-1)+1,nd*(lm(j)-1)+2];
        else
            Nf1=nd*N1;
            lmp=lm(j)-N1;
            lmg=[lmg,Nf1+2*lmp-1,Nf1+2*lmp];  %degrees of freedom indices
        end
    end
    
    for i=1:length(gpWts)    %loop over all Gauss points
        s = gpLocs(i, 1); h = gpLocs(i, 2); w = gpWts(i);
        n = [1-h-s,s,h];
        dnh=[-1,0,1];
        dns=[-1,1,0];
        
        dxs = dns*p(lm,1); dxh = dnh*p(lm,1);
        dys = dns*p(lm,2); dyh = dnh*p(lm,2);     
        J = [dxs, dys; dxh, dyh]; %pix
        detJ = det(J);  %pix^2
        dnx = (J(2, 2)*dns - J(1, 2)*dnh)/detJ;  %1/pix
        dny = (-J(2, 1)*dns + J(1, 1)*dnh)/detJ;

        b0=zeros(3,6);      
        
        if par.stiffnessgradient
            %local coordinates:
            xl=n*p(lm,1);
            yl=n*p(lm,2);
            dr=norm(par.rnotch-[xl,yl]);
            Ecor=1+(par.dE/par.E)*exp(-(dr/par.decaylength).^par.exp);
        else
            Ecor=1;
        end

        for j=1:3
            b0(1,2*j-1)=dnx(j);  %1/pix
            b0(2,2*j)=dny(j);
            b0(3,2*j-1)=dny(j);
            b0(3,2*j)=dnx(j);   
        end
        
        kl = kl + Ecor*0.5*detJ*par.h*w*b0'*C*b0;  %int b'*C*b dV  
    end
    for krow=1:6  % Here, we use a fast algorithm for sparse matrix assembly
        for kcol=1:6
            if abs(kl(krow,kcol))>0
                ntriplets = ntriplets + 1 ;
                len = length (KX) ;
                if (ntriplets > len)  % grow arrays
                    KI (2*len) = 0 ;
                    KJ (2*len) = 0 ;
                    KX (2*len) = 0 ;
                end
                KI (ntriplets) = lmg(krow);
                KJ (ntriplets) = lmg(kcol) ;
                KX (ntriplets) = kl(krow,kcol) ;
            end
        end
    end
    
end  %end loop over all elements
K= sparse(KI(1:ntriplets),KJ(1:ntriplets),KX(1:ntriplets),Nf,Nf);
