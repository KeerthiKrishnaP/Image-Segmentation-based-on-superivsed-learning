function saveVoltif(VOL0,nfich)
% save the 3d image into tif format of 8 bit
%   VOL0: the volume to save
%   nfich: the saving file name


disp(['-->Saving',nfich, 'in tif-file']);
nB = size(VOL0,3);
for i=1:nB
    imwrite(VOL0(:,:,i),nfich,'WriteMode','append');
end
disp(['-->Saving completed']);

end
