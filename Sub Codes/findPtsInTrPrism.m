function [id] = findPtsInTrPrism(xl,yl,zl,l1,l2,l3,vtx,n)
%find points within the triangular prism that is obtained by extrusion of a
%triangle
%
%
% Note: the three side vectors must be calculated by following the order:
%     l1 = vtx(1,:)-vtx(2,:);
%     l2 = vtx(2,:)-vtx(3,:);
%     l3 = vtx(3,:)-vtx(1,:);

% plane of the first side
nn = cross(l1,n);
id1 = nn(1).*(xl-vtx(1,1)) + ...
      nn(2).*(yl-vtx(1,2)) + ...
      nn(3).*(zl-vtx(1,3));
cr = nn(1).*(vtx(3,1)-vtx(1,1)) + ...
     nn(2).*(vtx(3,2)-vtx(1,2)) + ...
     nn(3).*(vtx(3,3)-vtx(1,3));
if cr>0
    id1 = id1 >= 0;
elseif cr<0
    id1 = id1 <= 0;
else
    error('Three vertices are on the same line !')
end

% plane of the second side
nn = cross(l2,n);
id2 = nn(1).*(xl-vtx(2,1)) + ...
      nn(2).*(yl-vtx(2,2)) + ...
      nn(3).*(zl-vtx(2,3));
cr = nn(1).*(vtx(1,1)-vtx(2,1)) + ...
     nn(2).*(vtx(1,2)-vtx(2,2)) + ...
     nn(3).*(vtx(1,3)-vtx(2,3));
if cr>0
    id2 = id2 >= 0;
elseif cr<0
    id2 = id2 <= 0;
else
    error('Three vertices are on the same line !')
end

% plane of the thrid side
nn = cross(l3,n);
id3 = nn(1).*(xl-vtx(3,1)) + ...
      nn(2).*(yl-vtx(3,2)) + ...
      nn(3).*(zl-vtx(3,3));
cr = nn(1).*(vtx(2,1)-vtx(3,1)) + ...
     nn(2).*(vtx(2,2)-vtx(3,2)) + ...
     nn(3).*(vtx(2,3)-vtx(3,3));
if cr>0
    id3 = id3 >= 0;
elseif cr<0
    id3 = id3 <= 0;
else
    error('Three vertices are on the same line !')
end

% 
id = id1 & id2 & id3;
