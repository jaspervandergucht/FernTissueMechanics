function [el,s,h] = FindEnclosingElement(t,p,coord)


%first find closest node:
distances = sum(bsxfun(@minus, p, coord).^2,2); %sqaured distance matrix
[closestdistancec,closestnode]=min(distances);

%find all elements connected to closestnode:
[ie,je]=find(t==closestnode);

el=0;
s=0;
h=0;
TOL=1e-6;
for i=1:numel(ie)
    e=ie(i);
    lm=t(e,:);  %node indices
    xy=p(lm,:);
    J=[xy(2,1)-xy(1,1),xy(3,1)-xy(1,1);xy(2,2)-xy(1,2),xy(3,2)-xy(1,2)];
    sh=J\[coord(1)-xy(1,1);coord(2)-xy(1,2)];
    if sh(1)>-TOL && sh(1)<1+TOL && sh(2)>-TOL && sh(2)<1+TOL && 1-sh(1)-sh(2)>-TOL && 1-sh(1)-sh(2)<1+TOL 
        s=sh(1);
        h=sh(2);
        el=e;
        break
    
    end
end