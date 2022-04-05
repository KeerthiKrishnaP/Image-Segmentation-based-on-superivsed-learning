function [vxlst] = rmVox3D(A,SLlst,e)
%choose the voxels to be removed by mouse click on each slice
%

s = 2;
es = e*s; %to enlarge the local volume around each 3D triangle

nz = length(SLlst);
siz = size(A);

% choose polygon vertices
x = [];
y = [];
z = [];
figure;
for i=SLlst
    SL = A(:,:,i);
    imshow(SL);title(['slice number: ',num2str(i)]);
    
    w=waitforbuttonpress;
    if w==1; continue; end %keypress -> next slice
    
    [xi,yi] = ginput;
    zi = zeros(size(xi)) + i;
    x = [x;xi];
    y = [y;yi];
    z = [z;zi];
end


% find the voxels within each triangle (distance < e)
% loop for each two adjacent slices
vxlst = cell(nz-1,1);
for i=1:nz-1
    % find the points at two adjacent slices
    id = z==SLlst(i);
    xi0 = x(id);
    yi0 = y(id);
    zi0 = z(id);
    id = z==SLlst(i+1);
    xi1 = x(id);
    yi1 = y(id);
    zi1 = z(id);
    
    np0 = length(xi0);
    np1 = length(xi1);
    if np0==0 || np1==0; continue; end
    
    % connectivity matrix of 3D triangles
    npmax = max(np0,np1);
    npmin = min(np0,np1);
    cr = (npmin-1)*2;
    conn = zeros(np0+np1-2,3,'uint32');
    ntr = size(conn,1);
    a = 1:npmin;
    b = npmin + [1:npmax];
    t = [npmax+1:npmax+npmin,1:npmax]; %potential use if np0>np1
    conn(1,:) = [a(1) b(1) b(2)];
    ai = 2;
    bi = 3;
    for j= 2:ntr
        if j<=cr && mod(j,2)==0
            conn(j,:) = [conn(j-1,[3,1]), a(ai)];
            ai = ai+1;
        elseif j<=cr && mod(j,2)~=0
            conn(j,:) = [conn(j-1,[3,1]), b(bi)];
            bi = bi+1;
        else
            conn(j,:) = [a(end),b(bi-1),b(bi)];
            bi = bi+1;
        end
    end
    if np0>np1 %if np0>np1, inverse the indices for a and b 
        conn = t(conn);
    end
    
    % loop for each 3D triangle
    pt = [[xi0;xi1],[yi0;yi1],[zi0;zi1]];
    vxlst0 = cell(ntr,1);
    for j = 1:ntr
        % coords of local volume around the triangle
        vtx = [pt(conn(j,:),:)];  %vertices of the triangle
        xmin = min(vtx(:,1))-es;  xmax = max(vtx(:,1))+es;
        ymin = min(vtx(:,2))-es;  ymax = max(vtx(:,2))+es;
        zmin = min(vtx(:,3))-es;  zmax = max(vtx(:,3))+es;
        if xmin<1; xmin=1; end;  if xmax>siz(2); xmax=siz(2); end
        if ymin<1; ymin=1; end;  if ymax>siz(1); ymax=siz(1); end
        if zmin<1; zmin=1; end;  if zmax>siz(3); zmax=siz(3); end
        [yl,xl,zl] = ndgrid(ymin:ymax,xmin:xmax,zmin:zmax);
        
        % the plane on which the triangle is lying
        l1 = vtx(1,:)-vtx(2,:); %side vectors of the triangle
        l2 = vtx(2,:)-vtx(3,:);
        l3 = vtx(3,:)-vtx(1,:);
        n = cross(l1,l2);         %normal vector of the plane
        
        % distance to the plane
        d = abs( n(1).*(xl-vtx(1,1)) + ...
                 n(2).*(yl-vtx(1,2)) + ...
                 n(3).*(zl-vtx(1,3)) ) ./ norm(n);

        % find the coords within the triangular prism
        id = findPtsInTrPrism(xl,yl,zl,l1,l2,l3,vtx,n);
        
        % find the coords within the 3D triangle (dist<e)
        id = id & d<=e;
        
        % record the coords to voxel list
        vxlst0{j} = sub2ind(siz,yl(id),xl(id),zl(id));
    end
    % record the coords to voxel list
    vxlst{i} = cell2mat(vxlst0);
end
vxlst = uint32(unique(cell2mat(vxlst)));

