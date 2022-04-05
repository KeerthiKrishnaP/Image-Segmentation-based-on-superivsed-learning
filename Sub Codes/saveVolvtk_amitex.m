function saveVolvtk_amitex(VOL0,nfich,type,choix,spacing)
% save the 3d image into vtk format
%
%       saveVolvtk(VOL0,nfich,type,choix)
%       saveVolvtk(VOL0,nfich,type,choix,spacing)
%       ---------------------------------
%   VOL0: the volume to save
%       note: in the case that more than one components of the volume
%             should be saved, VOL0 is cell data. 
%   nfich: the saving file name
%   type: 'logical', 'uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32',
%           'uint64', 'int64', 'single', 'double'
%   choix: 'point', 'cell' corresponding resp. {POINT_DATA}, {CELL_DATA}
%   spacing: (vector of 3x1) the spacing of regular grid
%

%  Keerthi Krishna PARVATHANENI  16/04/2018
%
% comment 16/04/2018 Yang Chen:
% the saveVolvtk_amitex is used for the application of AMITEX software,
% which follows the same conventionale X Y as MATLAB, so no need to permute
% the first two dimensions of the matrix volume

dx=1;dy=1;dz=1;
if nargin>=5
    dx=spacing(1); dy=spacing(2); dz=spacing(3);
end

tic
disp(['saving the volume in vtk-file: ',nfich]);

% ------ if input VOL0 is not cell data => only one component to be saved
if ~iscell(VOL0)
    disp 'input volume contains only one component...'
%     [nx, ny, nz] = size(VOL0);
    [nx, ny, nz] = size(VOL0);
    fid = fopen(nfich, 'w');
    fprintf(fid, '# vtk DataFile Version 4.5\n');
    fprintf(fid, 'tomo_Volume\n');
    fprintf(fid, 'BINARY\n');
    fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
    % CELL_DATA or POINT_DATA
    if strcmp(choix,'cell')
        fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', nx+1, ny+1, nz+1);
        fprintf(fid, 'ORIGIN    0.000   0.000   0.000\n');
        fprintf(fid, 'SPACING    %e    %e   %e\n', dx, dy, dz);
        fprintf(fid, 'CELL_DATA   %lu\n', nx*ny*nz);
    elseif strcmp(choix,'point')
        fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', nx, ny, nz);
        fprintf(fid, 'ORIGIN    0.000   0.000   0.000\n');
        fprintf(fid, 'SPACING    %f    %f   %f\n', dx, dy, dz);
        fprintf(fid, 'POINT_DATA   %lu\n', nx*ny*nz);
    end
    % dataType
    if (strcmp(type,'logical')==1)
        fprintf(fid, 'SCALARS tomo_Volume bit\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        fwrite(fid,VOL0(:),'ubit1',0,'b');
    elseif (strcmp(type,'uint8')==1)
        fprintf(fid, 'SCALARS tomo_Volume unsigned_char\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        fwrite(fid,VOL0(:),'uint8',0,'b');
    elseif (strcmp(type,'int8')==1)
        fprintf(fid, 'SCALARS tomo_Volume char\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        fwrite(fid,VOL0(:),'int8',0,'b');
    elseif (strcmp(type,'uint16')==1)
        fprintf(fid, 'SCALARS tomo_Volume unsigned_short\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        fwrite(fid,VOL0(:),'uint16',0,'b');
    elseif (strcmp(type,'int16')==1)
        fprintf(fid, 'SCALARS tomo_Volume short\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        fwrite(fid,VOL0(:),'int16',0,'b');
    elseif (strcmp(type,'uint32')==1)
        fprintf(fid, 'SCALARS tomo_Volume unsigned_int\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        fwrite(fid,VOL0(:),'uint32',0,'b');
    elseif (strcmp(type,'int32')==1)
        fprintf(fid, 'SCALARS MaterialId int\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        fwrite(fid,VOL0(:),'int32',0,'b');
    elseif (strcmp(type,'uint64')==1)
        fprintf(fid, 'SCALARS MaterialId unsigne_long\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        fwrite(fid,VOL0(:),'uint64',0,'b');
    elseif (strcmp(type,'int64')==1)
        fprintf(fid, 'SCALARS MaterialId long\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        fwrite(fid,VOL0(:),'int64',0,'b');
    elseif (strcmp(type,'single')==1)
        fprintf(fid, 'SCALARS MaterialId float\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        fwrite(fid,VOL0(:),'float',0,'b');
    elseif (strcmp(type,'double')==1)
        fprintf(fid, 'SCALARS MaterialId double\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        fwrite(fid,VOL0(:),'double',0,'b');
    end
    fclose(fid);
% ----------

% ------ if input VOL0 is cell data => several components to be saved
elseif iscell(VOL0)
    disp 'input volume contains several components...'
    % header of the file
%     [nx, ny, nz] = size(VOL0{1});
    [nx, ny, nz] = size(VOL0{1});
    fid = fopen(nfich, 'w');
    fprintf(fid, '# vtk DataFile Version 4.5\n');
    fprintf(fid, 'tomo_Volume\n');
    fprintf(fid, 'BINARY\n');
    fprintf(fid, 'DATASET STRUCTURED_POINTS\n');

    % CELL_DATA or POINT_DATA
    if strcmp(choix,'cell')
        fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', nx+1, ny+1, nz+1);
        fprintf(fid, 'ORIGIN    0.000   0.000   0.000\n');
        fprintf(fid, 'SPACING    %f    %f   %f\n', dx, dy, dz);
        fprintf(fid, 'CELL_DATA   %lu\n', nx*ny*nz);
    elseif strcmp(choix,'point')
        fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', nx, ny, nz);
        fprintf(fid, 'ORIGIN    0.000   0.000   0.000\n');
        fprintf(fid, 'SPACING    %f    %f   %f\n', dx, dy, dz);
        fprintf(fid, 'POINT_DATA   %lu\n', nx*ny*nz);
    end
    % dataType
    if (strcmp(type,'logical')==1)
        typeVTK = 'bit'; type = 'uint1';
    elseif (strcmp(type,'uint8')==1)
        typeVTK = 'unsigned_char';
    elseif (strcmp(type,'int8')==1)
        typeVTK = 'char';
    elseif (strcmp(type,'uint16')==1)
        typeVTK = 'unsigned_short';
    elseif (strcmp(type,'int16')==1)
        typeVTK = 'short';
    elseif (strcmp(type,'uint32')==1)
        typeVTK = 'unsigned_int';
    elseif (strcmp(type,'int32')==1)
        typeVTK = 'int';
    elseif (strcmp(type,'uint64')==1)
        typeVTK = 'unsigne_long';
    elseif (strcmp(type,'int64')==1)
        typeVTK = 'long';
    elseif (strcmp(type,'single')==1)
        typeVTK = 'float';
    elseif (strcmp(type,'double')==1)
        typeVTK = 'double';
    end
    % LOOKUP table
    nfields = length(VOL0);
    disp(['Number of fields to be saved: ',num2str(nfields)]);
    for ifield = 1:nfields
        fprintf(fid, ['SCALARS field_',num2str(ifield),' ',typeVTK,'\n']);
        fprintf(fid, 'LOOKUP_TABLE default\n');
        fwrite(fid,VOL0{ifield}(:),type,0,'b');
    end
    
    fclose(fid);
end
% ----
toc
end
