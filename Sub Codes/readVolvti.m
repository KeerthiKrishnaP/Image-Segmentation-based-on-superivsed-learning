function [ X, dx, dy, dz, typeMAT, name ] = readVolvti( varargin )
%Read a vti volume file
%
%   [ X ] = readVolvti( nfich )
%   [ X ] = readVolvti( nfich, nblocs )
%   [ X ] = readVolvti( nfich, nblocs, [n1, n2] )
% ========================================================
%
% Inputs:
%   - nfich: file name (*.vtk)
%   - nblocs: it is possible to read a huge volume by nblocs times =>
%       reduce the memory to be used
%   - [n1, n2]: the volume is read from n1-th slice to n2-th slice
% -----------------------------------------------------------------------
%
% Outputs:
%   - X: Volume data in the reading file (nfich)
% ====================================================================
%
% Keerthi Krishna PARVATHANENI  14/11/2016

% 
% ================= preparation
% Check the number of input arguments.
narginchk(1,3);
%
switch nargin   
    case 1
       nfich = varargin{1};
       nblocs=1;
    case 2
       nfich = varargin{1};
       nblocs = varargin{2};
    case 3
       nfich = varargin{1};
       nblocs = varargin{2};
       limits = varargin{3};
       n1 = limits(1);    n2 = limits(2);
end

tic
disp(['reading the volume from vti-file: ',nfich]);

% =======================================
fid=fopen(nfich,'r');

A='aaaaaaaaaa';
while (strcmp(A,'<ImageData')==0)
    ligne = fgetl(fid);
    ligne = strtrim(ligne);
    if (size(ligne,2)>=10);  A = ligne(1:10);  end
end

k1 = findstr(ligne,'WholeExtent');
k2 = findstr(ligne,'Origin');
k3 = findstr(ligne,'Spacing');

% DIMENSION
k11 = findstr(ligne(k1:k2),'"');
tmp = ligne(k1+k11(1):k1+k11(2)-2);
tmp = str2num(tmp);
ny = tmp(2);nx=tmp(4);nz=tmp(6);
if nargin==1 || nargin==2
    n1=1;  n2=nz;
end

% ORIGIN
k11 = findstr(ligne(k2:k3),'"');
tmp = ligne(k2+k11(1):k2+k11(2)-2);
tmp = str2num(tmp);
Oy = tmp(1);Ox=tmp(2);Oz=tmp(3);

% SPACING
k11 = findstr(ligne(k3:end),'"');
tmp = ligne(k3+k11(1):k3+k11(2)-2);
tmp = str2num(tmp);
dy = tmp(1);dx=tmp(2);dz=tmp(3);

disp(['sizeVOL = ',num2str([nx ny nz])]); 

A='aaaaaaaaa';
while (strcmp(A,'<CellData')==0)
    ligne = fgetl(fid);
    ligne = strtrim(ligne);
    if (size(ligne,2)>=9);  A = ligne(1:9);  end
end
ligne = fgetl(fid);  ligne = strtrim(ligne);
if ~strcmp(ligne,'</CellData>')
    disp('Cell data reading')
    dataType = 'CellData';
    nz = n2-n1+1;
    error('Cell data reading is not yet available !')
end

A='aaaaaaaaaa';
while (strcmp(A,'<PointData')==0)
    ligne = fgetl(fid);
    ligne = strtrim(ligne);
    if (size(ligne,2)>=10);  A = ligne(1:10);  end
end
ligne = fgetl(fid);  ligne = strtrim(ligne);
if ~strcmp(ligne,'</PointData>')
    disp('Point data reading')
    dataType = 'PointData';
    nx = nx+1;  ny = ny+1;  nz = n2-n1+1;
    k1 = findstr(ligne,'type');
    k2 = findstr(ligne,'Name');
    k3 = findstr(ligne,'format');
    k4 = findstr(ligne,'offset');
    
    % DATA FORMAT
    k11 = findstr(ligne(k1:k2),'"');
    typeVTI = ligne(k1+k11(1):k1+k11(2)-2);
    
    % NAME
    k11 = findstr(ligne(k2:k3),'"');
    name = ligne(k2+k11(1):k2+k11(2)-2);
    
    % WRITING FORMAT
    k11 = findstr(ligne(k3:k4),'"');
    format = ligne(k3+k11(1):k3+k11(2)-2);
    
    % OFFSET
    k11 = findstr(ligne(k4:end),'"');
    offset = ligne(k4+k11(1):k4+k11(2)-2);
    offset = str2num(offset);
end

A='aaaaaaaaaaaa';
while (strcmp(A,'</ImageData>')==0)
    ligne = fgetl(fid);
    ligne = strtrim(ligne);
    if (size(ligne,2)>=12);  A = ligne(1:12);  end
end

% ENCODING
ligne = fgetl(fid);    ligne = strtrim(ligne);
k1 = findstr(ligne,'encoding');
k11 = findstr(ligne(k1:end),'"');
encoding = ligne(k1+k11(1):k1+k11(2)-2);

% DATA READING

fgets(fid,3);% move to the the begining of the data
ntot = fread(fid,1,'uint64',0,'b');

typeMAT = typeVTI2MAT(typeVTI);

nbytes = (n1-1) * (nx*ny) * nBytes(typeMAT);
fseek(fid,nbytes,'cof');% move to the the begining slice (n1-th)

X = zeros(nx,ny,nz,typeMAT);

inc=ceil(nz/nblocs);  over=0;
for ibloc=1:nblocs
    i1 = 1+inc*(ibloc-1);
    i2 = i1+inc-1;
    if i2>nz; i2=nz; inc=i2-i1+1; over=1;  end
    B=fread(fid,nx*ny*inc,typeMAT,0,'b');
    X(:,:,i1:i2)=permute(reshape(B,ny,nx,inc),[2 1 3]);%"permute" to adjust the fread with the conventional coordinates (x<=>column, y<=>row)
    disp(['bloc',num2str(ibloc),': ',num2str([i1+n1-1,i2+n1-1])])
    if over==1; break; end
end

% if strcmp(dataType,'PointData')
%     X(end,:,:) = [];
%     X(:,end,:) = [];
%     X(:,:,end) = [];
% end

fclose(fid);

toc


%% functions
function [typeMAT] = typeVTI2MAT(typeVTI)
if (strcmp(typeVTI,'UInt8')==1)
    typeMAT = 'uint8';
elseif (strcmp(typeVTI,'Int8')==1)
    typeMAT = 'int8';
elseif (strcmp(typeVTI,'UInt16')==1)
    typeMAT = 'uint16';
elseif (strcmp(typeVTI,'Int16')==1)
    typeMAT = 'int16';
elseif (strcmp(typeVTI,'UInt32')==1)
    typeMAT = 'uint32';
elseif (strcmp(typeVTI,'Int32')==1)
    typeMAT = 'int32';
elseif (strcmp(typeVTI,'UInt64')==1)
    typeMAT = 'uint64';
elseif (strcmp(typeVTI,'Int64')==1)
    typeMAT = 'int64';
elseif (strcmp(typeVTI,'Float32')==1)
    typeMAT = 'single';
elseif (strcmp(typeVTI,'Float64')==1)
    typeMAT = 'double';
end
end


function [nbytes] = nBytes(typeMAT)
if (strcmp(typeMAT,'uint8')==1)
    nbytes = 1;
elseif (strcmp(typeMAT,'int8')==1)
    nbytes = 1;
elseif (strcmp(typeMAT,'uint16')==1)
    nbytes = 2;
elseif (strcmp(typeMAT,'int16')==1)
    nbytes = 2;
elseif (strcmp(typeMAT,'uint32')==1)
    nbytes = 4;
elseif (strcmp(typeMAT,'int32')==1)
    nbytes = 4;
elseif (strcmp(typeMAT,'uint64')==1)
    nbytes = 8;
elseif (strcmp(typeMAT,'int64')==1)
    nbytes = 8;
elseif (strcmp(typeMAT,'single')==1)
    nbytes = 4;
elseif (strcmp(typeMAT,'float')==1)
    nbytes = 4;
elseif (strcmp(typeMAT,'double')==1)
    nbytes = 8;
end
end



end
