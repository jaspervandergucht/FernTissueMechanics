function [R] = BoundaryStresses(p,eb,par,Nf)
%Obtain nodal forces and moments due to the pressure par.P at the boundary

if nargin<4
    if par.includebending
        Nf=3*size(p,1);
    else
        Nf=2*size(p,1);
    end
end

R=zeros(Nf,1);
for i=1:size(eb,1);
    lm=eb(i,:);
    x1=p(lm(1),1);y1=p(lm(1),2);
    x2=p(lm(2),1);y2=p(lm(2),2);
    L=sqrt((x2-x1)^2+(y2-y1)^2);
    c=(x2-x1)/L; %direction cosine
    s=(y2-y1)/L; %direction sine

    if par.includebending
        lmg=[3*(lm(1)-1)+1,3*(lm(1)-1)+2,3*(lm(1)-1)+3,3*(lm(2)-1)+1,3*(lm(2)-1)+2,3*(lm(2)-1)+3];
        ri=0.5*par.P*L*[0;1;L/6;0;1;-L/6];
        T=[c,s,0,0,0,0;
           -s,c,0,0,0,0;
            0,0,1,0,0,0;
            0,0,0,c,s,0;
            0,0,0,-s,c,0;
            0,0,0,0,0,1];
        R(lmg)=R(lmg)+T'*ri;
    else
        lmg=[2*(lm(1)-1)+1,2*(lm(1)-1)+2,2*(lm(2)-1)+1,2*(lm(2)-1)+2];
        ri=0.5*par.P*L*[0;1;0;1];
        T=[c,s,0,0;
           -s,c,0,0;
            0,0,c,s;
            0,0,-s,c];
        R(lmg)=R(lmg)+T'*ri;
    end

    
    
    
end