%--------------------------------------------------------------------------
%
% FONCTION DE SAUVEGARDE AU FORMAT VTI, en POINTDATA
% d'un champ d'entier (défini sur la base de ndgrid)
%
% function sauv_vti(X,dx,dy,dz,nomfic,type,S)
%
%       X : 		champ a sauvegarder (issu de ndgrid)
%	dx,dy,dz : 	dimensions d'un voxel (spacing)
%	nomfic :	nom de fichier 
%	S :		nom du champ
%
%
% FORMAT VTI : PointData, type, Big Endian
%	type = 'string' (type VTK, type VTI)
%	type = 'uint8' (unsigned_char,UInt8)
%	type = 'uint16' (unsigned_short,UInt16)
%	type = 'uint32' (unsigned_int,UInt32)
%	type = 'uint64' (unsigned_long,UInt64)
%	type = 'int8' (char, Int8)
%	type = 'int16' (short,Int16)
%	type = 'int32' (int,Int32)
%	type = 'int64' (long,Int64)
%       type = 'float' (float,Float32)
%
%--------------------------------------------------------------------------
function sauv_vti(X,dx,dy,dz,nomfic,type,S)

if (nargin~=7);error('nombre d arguments incorrect');end;

if (strcmp(type,'uint8')==1);
  k=1;
  typeVTK='UInt8';
elseif (strcmp(type,'uint16')==1);
  k=2;
  typeVTK='UInt16';
elseif (strcmp(type,'uint32')==1);
  k=4;
  typeVTK='UInt32';
elseif (strcmp(type,'uint64')==1);
  k=8;
  typeVTK='UInt64';
elseif (strcmp(type,'int8')==1);
  k=1;
  typeVTK='Int8';
elseif (strcmp(type,'int16')==1);
  k=2;
  typeVTK='Int16';
elseif (strcmp(type,'int32')==1);
  k=4;
  typeVTK='Int32';
elseif (strcmp(type,'int64')==1);
  k=8;
  typeVTK='Int64';
elseif (strcmp(type,'float')==1);
  k=4;
  typeVTK='Float32';
else
error('type invalide');
end;

disp(['saving vti-file: ',nomfic])

X = permute(X,[2,1,3]);%"permute" to adjust the fread with the conventional coordinates (x<=>column, y<=>row)

%----ECRITURE DEBUT DU FICHIER VTI : adaptation pointdata
    [nx, ny, nz] = size(X);
    fid = fopen(nomfic, 'w');
    fprintf(fid, '<?xml version="1.0"?>\n');
    fprintf(fid, '<VTKFile type="ImageData" version="1.0" byte_order="BigEndian" header_type="UInt64">\n');
%    fprintf(fid, '  <ImageData WholeExtent="0 %d 0 %d 0 %d" Origin="0 0 0" Spacing="%f %f %f">\n',nx,ny,nz,dx,dy,dz);
    fprintf(fid, '  <ImageData WholeExtent="0 %d 0 %d 0 %d" Origin="0 0 0" Spacing="%f %f %f">\n',nx-1,ny-1,nz-1,dx,dy,dz);
    fprintf(fid, '    <Piece Extent="0 %d 0 %d 0 %d">\n',nx-1,ny-1,nz-1);
%    fprintf(fid, '      <PointData>\n');
%    fprintf(fid, '      </PointData>\n');
    fprintf(fid, '      <CellData>\n');
    fprintf(fid, '      </CellData>\n');
%    fprintf(fid, '      <CellData Scalars="%s">\n',S);
    fprintf(fid, '      <PointData Scalars="%s">\n',S);
    fprintf(fid, '        <DataArray type="%s" Name="%s" format="appended" offset="0"/>\n',typeVTK,S);
%    fprintf(fid, '      </CellData>\n');
    fprintf(fid, '      </PointData>\n');
    fprintf(fid, '    </Piece>\n');
    fprintf(fid, '  </ImageData>\n');
    fprintf(fid, '  <AppendedData encoding="raw">\n');
    fprintf(fid, '  _');

%-----ECRITURE DES DONNEES
    fwrite(fid,k*nx*ny*nz,'uint64',0,'b');
    fwrite(fid,X,type,0,'b');

%----FIN DES DONNEES
    fprintf(fid, '\n  </AppendedData>\n');
    fprintf(fid, '</VTKFile>\n');

fclose(fid);
return;

