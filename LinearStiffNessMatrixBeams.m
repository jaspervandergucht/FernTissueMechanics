function [K] = LinearStiffNessMatrixBeams(p,e,par,Nf) 
    
% Calculate stiffness matrix, considering each element in e as a 2-node spring with
% stiffness EA and bending rigidity EI
% p: node coordinates
% e: edge-list
% par contains E( Young modulus), A (beam cross section), I (moment of inertia) 
% Nf: # of degrees of freedom (if not specified, equal to 3*size(p,1)
% Output: 
% K: stiffness matrix


% Use sparse matrix utility, see:
% http://blogs.mathworks.com/loren/2007/03/01/creating-sparse-finit
% e-element-matrices-in-matlab/
if nargin<4
    Nf=3*size(p,1);
end
I = zeros (10, 1) ;
J = zeros (10, 1) ;
X = zeros (10, 1) ;
ntriplets = 0 ;
 
% loop over elements in e:
for i=1:size(e,1)
    lm=e(i,:);
    lmg=[3*(lm(1)-1)+1,3*(lm(1)-1)+2,3*(lm(1)-1)+3,3*(lm(2)-1)+1,3*(lm(2)-1)+2,3*(lm(2)-1)+3];

    x1=p(lm(1),1);y1=p(lm(1),2);
    x2=p(lm(2),1);y2=p(lm(2),2);
    L=sqrt((x2-x1)^2+(y2-y1)^2);
    c=(x2-x1)/L; %direction cosine
    s=(y2-y1)/L; %direction sine
    if par.stiffnessgradient
        rmid=0.5*[x1+x2,y1+y2];
        dr=norm(rmid-par.rnotch);
        E=1+(par.dE/par.E)*exp(-(dr/par.decaylength).^par.exp)
        
    else
        E=par.E;
    end
    %L=L*par.pixelsize;
    ks=E*par.A/(L);
    kb=E*par.I/L^3;

    ki=[ks, 0, 0, -ks , 0, 0;
        0, 12*kb, 6*kb*L, 0, -12*kb, 6*kb*L;
        0, 6*kb*L, 4*kb*L^2, 0, -6*kb*L,2*kb*L^2;
        -ks, 0,0, ks , 0 ,0;
        0, -12*kb, -6*kb*L,0,12*kb,-6*kb*L;
        0, 6*kb*L, 2*kb*L^2,0,-6*kb*L,4*kb*L^2];


     % (Book bhatti, p 266 )  
    T=[c,s,0,0,0,0;
       -s,c,0,0,0,0;
        0,0,1,0,0,0;
        0,0,0,c,s,0;
        0,0,0,-s,c,0;
        0,0,0,0,0,1];
     
    kig=T'*ki*T;   %element stiffness in global coordinates
    
    
    for krow=1:6  % Here, we use a fast algorithm for sparse matrix assembl
        for kcol=1:6
            ntriplets = ntriplets + 1 ;
            len = length (X) ;
            if (ntriplets > len)  % grow arrays
                I (2*len) = 0 ;
                J (2*len) = 0 ;
                X (2*len) = 0 ;
            end
                I (ntriplets) = lmg(krow);
                J (ntriplets) = lmg(kcol) ;
                X (ntriplets) = kig (krow,kcol) ;
        end
    end
    


end

K = sparse (I(1:ntriplets),J(1:ntriplets),X(1:ntriplets),Nf,Nf) ;