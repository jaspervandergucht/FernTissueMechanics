function [K] = LinearStiffNessMatrixSprings(p,e,par,Nf) 
    
% Calculate stiffness matrix, considering each element in e as a 2-node spring with
% stiffness EA
% p: node coordinates
% e: edge-list
% par contains E( Young modulus), A (beam cross section)
% Nf: # of degrees of freedom (if not specified, equal to 2*size(p,1)

% Output: 
% K: stiffness matrix as sparse matrix

% Use sparse matrix utility, see:
% http://blogs.mathworks.com/loren/2007/03/01/creating-sparse-finit
% e-element-matrices-in-matlab/

if nargin<4
    Nf=2*size(p,1);
end

I = zeros (10, 1) ;
J = zeros (10, 1) ;
X = zeros (10, 1) ;
ntriplets = 0 ;
 


% loop over elements in e:
for i=1:size(e,1)
    lm=e(i,:);
    lmg=[2*(lm(1)-1)+1,2*(lm(1)-1)+2,2*(lm(2)-1)+1,2*(lm(2)-1)+2];

    x1=p(lm(1),1);y1=p(lm(1),2);
    x2=p(lm(2),1);y2=p(lm(2),2);
    L=sqrt((x2-x1)^2+(y2-y1)^2);
    c=(x2-x1)/L; %direction cosine
    s=(y2-y1)/L; %direction sine
    
    if par.stiffnessgradient
        rmid=0.5*[x1+x2,y1+y2];
        dr=norm(rmid-par.rnotch);
        E=par.E+par.dE*exp(-(dr/par.decaylength).^par.exp);
    else
        E=par.E;
    end

    ks=E*par.A/(L);   %N/micron
    
    ki=ks*[c^2, c*s, -c^2, -c*s;
           c*s, s^2, -c*s, -s^2;
           -c^2, -c*s, c^2, c*s;
           -c*s, -s^2, c*s, s^2];
    for krow=1:4  % Here, we use a fast algorithm for sparse matrix assembl
        for kcol=1:4
            ntriplets = ntriplets + 1 ;
            len = length (X) ;
            if (ntriplets > len)  % grow arrays
                I (2*len) = 0 ;
                J (2*len) = 0 ;
                X (2*len) = 0 ;
            end
                I (ntriplets) = lmg(krow);
                J (ntriplets) = lmg(kcol) ;
                X (ntriplets) = ki (krow,kcol) ;
        end
    end
    


end

K = sparse (I(1:ntriplets),J(1:ntriplets),X(1:ntriplets),Nf,Nf) ;