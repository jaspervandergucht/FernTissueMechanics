function [results] = PostProcess(p,p0,theta,t,e,par)

%calculate stresses, energies etc:
% elementstress: stress tensor [sxx,syy,sxy] in each element
% principalstresses: principal stresses smax,smin in each elememt
% principaldirection: principal direction vector in each element
% nodalstresses; same but in nodes
% nodalprincipalstresses
%
% bendingenergy: bending energy for each wall segment
% stretchingenergy: stretch energy for each segment
% relativeelongation: (L-L0)/L0 for each wall segment

%Post-processing
Ne=size(e,1);
Ebend=zeros(Ne,1);
Estretch=zeros(Ne,1);
Elongation=zeros(Ne,1);

ks=par.E*par.A;
EI=par.E*par.I;
for i=1:Ne
    lm=e(i,:);
    
    x1=p(lm(1),1);y1=p(lm(1),2);
    x2=p(lm(2),1);y2=p(lm(2),2);
    L=sqrt((x2-x1)^2+(y2-y1)^2);

    x10=p0(lm(1),1);y10=p0(lm(1),2);
    x20=p0(lm(2),1);y20=p0(lm(2),2);
    L0=sqrt((x20-x10)^2+(y20-y10)^2);
    
    

    Elongation(i)=(L-L0)/L0;
    Estretch(i)=0.5*ks*(L-L0)^2;
    
    ui=p(lm,:)'-p0(lm,:)';
    ui=ui(:);

    c=(x2-x1)/L; %direction cosine
    s=(y2-y1)/L; %direction sine

    T=[c,s,0,0;
        -s,c,0,0;
        0,0,c,s;
        0,0,-s,c];
    vi=T*ui;
    v1=vi(2);
    v2=vi(4);
    if par.includebending
    Q1=theta(lm(1));
    Q2=theta(lm(2));
    Ebend(i)=EI*((2*Q1^2+2*Q2^2+Q1*Q2)/L0+6*(Q1+Q2)*(v1-v2)/L0^2+6*(v1-v2)^2/L0^3);
    end
end



switch (par.type)  %1. plane stress; 2. plane strain
case 1
    C = par.E/(1 - par.v^2)*[1, par.v, 0; par.v, 1, 0; 0, 0, (1 - par.v)/2];
case 2
    C = par.E/((1 + par.v)*(1 - 2*par.v))*[1 - par.v, par.v, 0; par.v, 1 - par.v, 0;
        0, 0, (1 - 2*par.v)/2];
end  

Nt=size(t,1);
S=zeros(Nt,3);%stress
P=zeros(Nt,2);
St=zeros(Nt,3);%strain
Pst=zeros(Nt,2);

u=zeros(2*size(p,1),1);
u(1:2:end)=p(:,1)-p0(:,1);
u(2:2:end)=p(:,2)-p0(:,2);

NS=zeros(size(p,1),3);
NSt=zeros(size(p,1),3);
NN=zeros(size(p,1),1);

gpLocs = [1/3,1/3];  
% Loop over all elements
for i=1:Nt
    lm=t(i,:);  %node indices
    lmg=[];
    for j=1:3
        lmg=[lmg,2*lm(j)-1,2*lm(j)];  %degrees of freedom indices
    end
    s = gpLocs(1); h = gpLocs(2); 
    n = [1-h-s,s,h];
    dnh=[-1,0,1];
    dns=[-1,1,0];
    
    dxs = dns*p(lm,1); dxh = dnh*p(lm,1);
    dys = dns*p(lm,2); dyh = dnh*p(lm,2);
    J = [dxs, dys; dxh, dyh]; detJ = det(J);
    dnx = (J(2, 2)*dns - J(1, 2)*dnh)/detJ;
    dny = (-J(2, 1)*dns + J(1, 1)*dnh)/detJ;

    b0=zeros(3,6);     

    for j=1:3
        b0(1,2*j-1)=dnx(j);
        b0(2,2*j)=dny(j);
        b0(3,2*j-1)=dny(j);
        b0(3,2*j)=dnx(j);   
    end

    if par.stiffnessgradient
            %local coordinates:
            xl=n*p(lm,1);
            yl=n*p(lm,2);
            dr=norm(par.rnotch-[xl,yl]);
            Ecor=1+(par.dE/par.E)*exp(-dr/par.decaylength);
        else
            Ecor=1;
        end


    Sti=b0*u(lmg);  %strain
    Si=Ecor*C*b0*u(lmg);

    S(i,:)=Si';
    St(i,:)=Sti';
    NS(lm,:)=NS(lm,:)+repmat(Si',3,1);
    NSt(lm,:)=NSt(lm,:)+repmat(Sti',3,1);
    NN(lm)=NN(lm)+1;
    
end  %end loop over all elements
NS=NS./repmat(NN,1,3);
NSt=NSt./repmat(NN,1,3);

%Principal stresses:
wortel=sqrt((S(:,1)-S(:,2)).^2+4*S(:,3).^2);
P=[(S(:,1)+S(:,2)+wortel)/2,(S(:,1)+S(:,2)-wortel)/2];

wortel2=sqrt((NS(:,1)-NS(:,2)).^2+4*NS(:,3).^2);
NP=[(NS(:,1)+NS(:,2)+wortel2)/2,(NS(:,1)+NS(:,2)-wortel2)/2];

V1=[-0.5*(-S(:,1)+S(:,2)-wortel),S(:,3)]; %direction of maximum tensile stress
nV1=sqrt(V1(:,1).^2+V1(:,2).^2);
V1(:,1)=V1(:,1)./nV1;
V1(:,2)=V1(:,2)./nV1;

%Principal strains:
wortel=sqrt((St(:,1)-St(:,2)).^2+4*St(:,3).^2);
Pst=[(St(:,1)+St(:,2)+wortel)/2,(St(:,1)+St(:,2)-wortel)/2];

wortel2=sqrt((NSt(:,1)-NSt(:,2)).^2+4*NSt(:,3).^2);
NPst=[(NSt(:,1)+NSt(:,2)+wortel2)/2,(NSt(:,1)+NSt(:,2)-wortel2)/2];

V1s=[-0.5*(-St(:,1)+St(:,2)-wortel),St(:,3)]; %direction of maximum tensile strain
nV1s=sqrt(V1s(:,1).^2+V1s(:,2).^2);
V1s(:,1)=V1s(:,1)./nV1s;
V1s(:,2)=V1s(:,2)./nV1s;


results.elementstress=S;
results.principalstresses=P;
results.principaldirection=V1;

results.elementstrain=St;
results.principalstrains=Pst;
results.principalstraindirection=V1s;

results.bendingenergy=Ebend;
results.stretchingenergy=Estretch;
results.relativeelongation=Elongation;
results.nodalstresses=NS;
results.nodalprincipalstresses=NP;

results.nodalstrains=NSt;
results.nodalprincipalstrains=NPst;
%