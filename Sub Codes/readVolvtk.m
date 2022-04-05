function [ X,spacing,fieldnames,typeMAT ] = readVolvtk( varargin )
%Read a vtk volume file
%
%   [ X ] = readVolvtk( nfich )
%   [ X ] = readVolvtk( nfich, nblocs )
%   [ X ] = readVolvtk( nfich, n, n_tot )
%   [ X ] = readVolvtk( nfich, n, n_tot, nblocs )
%   [ X ] = readVolvtk( nfich, n, n_tot, nblocs, [n1, n2] )
%   [ X,spacing ] = readVolvtk( ... )
%   [ X,spacing, fieldname ] = readVolvtk( ... )
%   [ X,spacing, fieldname, typeMAT ] = readVolvtk( ... )
% ========================================================
%
% Inputs:
%   - nfich: file name (*.vtk)
%
%   - n: (vector) number of the components to be read in the volume file
%       example: n=[1 3 6] means to read the 1st, 3d and 6d components
%
%   - n_tot: total number of compnents contained in the volume "nfich"
%
%   - nblocs: it is possible to read a huge volume by nblocs times =>
%       reduce the memory to be used
%   - [n1, n2]: the volume is read from n1-th slice to n2-th slice
% -----------------------------------------------------------------------
%
% Outputs:
%   - X: Volume data in the reading file (nfich)
%       > if nargin==1 or 2, X is a array
%       > if nargin==3 or 4, X is a structure
%   - spacing:(vector3x1) the dimension spacing of the volume to be read
%   -fieldname: the name of the data to be read
%       > if nargin==1 or 2, fieldname is a string
%       > if nargin==3 or 4, fieldname is a cell of strings
%   - typeMAT: data types of fields in the volume to be read
%       > if nargin==1 or 2, typeMAT is a string
%       > if nargin==3 or 4, typeMAT is a cell of strings
% ====================================================================
%
% Keerthi Krishna PARVATHANENI  21/03/2016 

% 
% ================= preparation
% Check the number of input arguments.
narginchk(1,5);

%
switch nargin   
    case 1
       nfich = varargin{1};
       n=1; n_tot=1;
       nblocs=1;
    case 2
       nfich = varargin{1};
       n=1; n_tot=1;
       nblocs = varargin{2};
    case 3
       nfich = varargin{1};
       n = varargin{2};
       n_tot = varargin{3};
       nblocs=1;
    case 4
       nfich = varargin{1};
       n = varargin{2};
       n_tot = varargin{3};
       nblocs = varargin{4};
    case 5
       nfich = varargin{1};
       n = varargin{2};
       n_tot = varargin{3};
       nblocs = varargin{4};
       if isempty(nblocs); nblocs=1; end
       limits = varargin{5};
       n1 = limits(1);    n2 = limits(2);
end

%tic
disp(['reading the volume from vtk-file: ',nfich]);
%disp '------------'
%disp 'If you want to read several components, you must be sure about their number in your volume to be read !'
%disp 'Otherwise, there will be no end of this lecture !!!!'
%disp ' '
% =======================================
fid=fopen(nfich,'r');

% DIMENSION
A='aaaaaaaaaa';
while (strcmp(A,'DIMENSIONS')==0)
    ligne = fgetl(fid);
    if (size(ligne,2)>=10);  A = ligne(1:10);  end
end
tmp=[];
tmp = sscanf(ligne(11:size(ligne,2)),'%d');
% nx = tmp(1);ny=tmp(2);nz=tmp(3);
ny = tmp(1);nx=tmp(2);nz=tmp(3);
disp(['sizeVOL = ',num2str([nx ny nz])]); 

% SPACING
A='aaaaaaa';
while (strcmp(A,'SPACING')==0)
    ligne = fgetl(fid);
    if (size(ligne,2)>=7);  A = ligne(1:7);  end
end
tmp=[];
tmp = sscanf(ligne(8:size(ligne,2)),'%f');
dx = tmp(1);dy=tmp(2);dz=tmp(3);
spacing=[dx,dy,dz];
% dataType
A = 'aaaaaaaaa';
while (strcmp(A,'POINT_DAT')==0 && strcmp(A,'CELL_DATA')==0)
    ligne = fgetl(fid);
    if (size(ligne,2)>=9); A = ligne(1:9);  end
end
tmp = [];
tmp = sscanf(ligne(1:4),'%s');
if strcmp(tmp,'CELL');  nx=nx-1; ny=ny-1; nz=nz-1; end

% read field data
% ------------------
if nargin==1 || nargin==2   % to be compatible to the utilization of the old version of this fucntion
    while (strcmp(A,'SCALARS')==0)
        ligne = fgetl(fid);
        if (size(ligne,2)>=7);  A = ligne(1:7);  end
    end;
    tmp=[];
    tmp=textscan(ligne(8:end),'%s');
    fieldnames = char(tmp{1}(1));
    typeVTK = char(tmp{1}(2));
    % LOOKUP table
    A='aaaaaa';
    while (strcmp(A,'LOOKUP')==0)
        ligne = fgetl(fid);
        if (size(ligne,2)>=6);  A = ligne(1:6);  end
    end
    typeMAT = typeVTK2MAT(typeVTK);
    if(strcmp(typeVTK,'ubit1'))
        X = false(nx,ny,nz);  % cause "zeros(*,*,*,'ubit1')" not work
    else
        X = zeros(nx,ny,nz,typeMAT);
    end
    inc=ceil(nz/nblocs); over=0;
    for ibloc=1:nblocs
        i1 = 1+inc*(ibloc-1);
        i2 = i1+inc-1;
        if i2>nz; i2=nz; inc=i2-i1+1; over=1;  end
        B=fread(fid,nx*ny*inc,typeMAT,0,'b');
        X(:,:,i1:i2)=permute(reshape(B,ny,nx,inc),[2 1 3]);%"permute" to adjust the fread with the conventional coordinates (x<=>column, y<=>row)
        %X(:,:,i1:i2)= reshape(B,ny,nx,inc);
        disp(num2str([ibloc,i1,i2,inc]))
        if over==1; break; end
    end
end

if nargin==3 || nargin==4
    nb_chmps = length(n);
    X = cell(nb_chmps,1);
    fieldnames = cell(nb_chmps,1);
    typeMAT = cell(nb_chmps,1);
    i_chmp = 1;
    for i=1:n_tot
        while (strcmp(A,'SCALARS')==0)
            ligne = fgetl(fid);
            if (size(ligne,2)>=7);  A = ligne(1:7);  end
        end;
        tmp=[];
        tmp=textscan(ligne(8:end),'%s');
        fieldnames{i_chmp} = char(tmp{1}(1));
        typeVTK = char(tmp{1}(2));
        % LOOKUP table
        A='aaaaaa';
        while (strcmp(A,'LOOKUP')==0)
            ligne = fgetl(fid);
            if (size(ligne,2)>=6);  A = ligne(1:6);  end
        end
        if i==n(i_chmp)
            disp(['reading the field: ',fieldnames{i_chmp},' ',typeVTK]);
           % -----------------------------------
            typeMAT{i_chmp} = typeVTK2MAT(typeVTK);
            if(strcmp(typeVTK,'ubit1'))
                X{i_chmp} = false(nx,ny,nz);  % cause "zeros(*,*,*,'ubit1')" not work
            else
                X{i_chmp} = zeros(nx,ny,nz,typeMAT{i_chmp});
            end
            % 
            inc=ceil(nz/nblocs);  over=0;
            for ibloc=1:nblocs
                i1 = 1+inc*(ibloc-1);
                i2 = i1+inc-1;
                if i2>nz; i2=nz; inc=i2-i1+1; over=1;  end
                B=fread(fid,nx*ny*inc,typeMAT{i_chmp},0,'b');
                X{i_chmp}(:,:,i1:i2)=permute(reshape(B,ny,nx,inc),[2 1 3]);%"permute" to adjust the fread with the conventional coordinates (x<=>column, y<=>row)
                %X{i_chmp}(:,:,i1:i2)= reshape(B,ny,nx,inc);
                disp(num2str([ibloc,i1,i2,inc]))
                if over==1; break; end
            end
            i_chmp = i_chmp + 1; if(i_chmp>nb_chmps); i_chmp=nb_chmps; end
            % ------------------------------------------------------
        end
    end
    %
    % regroup the fields into structure data
    X = cell2struct(X,fieldnames,1);
end



if nargin==5
    nz = n2-n1+1;
    
    nb_chmps = length(n);
    X = cell(nb_chmps,1);
    fieldnames = cell(nb_chmps,1);
    typeMAT = cell(nb_chmps,1);
    i_chmp = 1;
    for i=1:n_tot
        while (strcmp(A,'SCALARS')==0)
            ligne = fgetl(fid);
            if (size(ligne,2)>=7);  A = ligne(1:7);  end
        end;
        tmp=[];
        tmp=textscan(ligne(8:end),'%s');
        fieldnames{i_chmp} = char(tmp{1}(1));
        typeVTK = char(tmp{1}(2));
        % LOOKUP table
        A='aaaaaa';
        while (strcmp(A,'LOOKUP')==0)
            ligne = fgetl(fid);
            if (size(ligne,2)>=6);  A = ligne(1:6);  end
        end
        if i==n(i_chmp)
            disp(['reading the field: ',fieldnames{i_chmp},' ',typeVTK]);
           % -----------------------------------
            typeMAT{i_chmp} = typeVTK2MAT(typeVTK);
            if(strcmp(typeVTK,'ubit1'))
                X{i_chmp} = false(nx,ny,nz);
            else
                X{i_chmp} = zeros(nx,ny,nz,typeMAT{i_chmp});
            end
            % 
            % move to the the begining slice (n1-th)
            nbytes = (n1-1) * (nx*ny) * nBytes(typeMAT{i_chmp});

            fseek(fid,nbytes,'cof');

            inc=ceil(nz/nblocs);  over=0;
            for ibloc=1:nblocs
                i1 = 1+inc*(ibloc-1);
                i2 = i1+inc-1;
                if i2>nz; i2=nz; inc=i2-i1+1; over=1;  end
                B=fread(fid,nx*ny*inc,typeMAT{i_chmp},0,'b');
                X{i_chmp}(:,:,i1:i2)=permute(reshape(B,ny,nx,inc),[2 1 3]);%"permute" to adjust the fread with the conventional coordinates (x<=>column, y<=>row)
                %X{i_chmp}(:,:,i1:i2)= reshape(B,ny,nx,inc);
                disp(num2str([ibloc,i1+n1-1,i2+n1-1,inc]))
                if over==1; break; end
            end
            i_chmp = i_chmp + 1; if(i_chmp>nb_chmps); i_chmp=nb_chmps; end
            % ------------------------------------------------------
        end
    end
    %
    % regroup the fields into structure data
    X = cell2struct(X,fieldnames,1);
end

fclose(fid);

%toc


function [typeMAT] = typeVTK2MAT(typeVTK)
if (strcmp(typeVTK,'bit')==1)
    typeMAT = 'ubit1';
elseif (strcmp(typeVTK,'unsigned_char')==1)
    typeMAT = 'uint8';
elseif (strcmp(typeVTK,'char')==1)
    typeMAT = 'int8';
elseif (strcmp(typeVTK,'unsigned_short')==1)
    typeMAT = 'uint16';
elseif (strcmp(typeVTK,'short')==1)
    typeMAT = 'int16';
elseif (strcmp(typeVTK,'unsigned_int')==1)
    typeMAT = 'uint32';
elseif (strcmp(typeVTK,'int')==1)
    typeMAT = 'int32';
elseif (strcmp(typeVTK,'unsigne_long')==1)
    typeMAT = 'uint64';
elseif (strcmp(typeVTK,'long')==1)
    typeMAT = 'int64';
elseif (strcmp(typeVTK,'float')==1)
    typeMAT = 'single';
elseif (strcmp(typeVTK,'double')==1)
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
