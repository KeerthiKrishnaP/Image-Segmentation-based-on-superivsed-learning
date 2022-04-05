function [varargout] = featureParam_bloc(A,r,str_out,nblocs)
%compute the feature parameters of a gray-level image by bloc
%in order to avoid the high usage of memory
%
%   [varargout] = featureParam_bloc(A,str_out,r,nblocs)
%   ---------------------------------------------------
%
%   Inputs:
%       > A : gray-level image to be processed
%       > str_out: string cell specifying the feature parameters to be
%                  computed (see help featureParam.m)
%       > r : parameter controlling the window size (cube of 2*r+1 / side)
%       > nblocs: number of blocs that the image will be cut into along its
%                 max dimension
%
%   Outputs:
%       > see help featureParam.m
%
% Keerthi Krishna PARVATHANENI 2018.04.04
%
% Abkup=A;
% A=Asub;

siz = size(A);
np = length(str_out);


% loop for bloc image
[dimMAX,dimId] = max(siz);
inc = ceil(dimMAX/nblocs); over=0;
for ibloc=1:nblocs
    % display
    disp(['--> bloc ',num2str(ibloc),' / ',num2str(nblocs)]);
    
    % cut the image into blocs along the max dimension
    i0 = 1 + (ibloc-1)*inc;
    i1 = i0 + inc - 1;
    if i1>dimMAX; i1=dimMAX; inc=i1-i0+1; over=1;  end
    
    % take some margins to avoid edge effects
    i01 = i0 - r-2;
    i11 = i1 + r+2;
    if i01<1; i01=1;end
    if i11>dimMAX; i11=dimMAX; end
    
    m0 = i0 - i01 + 1;
    m1 = i11 - i1;
    
    % find the current image bloc and its index in the whole image
    if dimId==1
        A_bloc = A(i01:i11,:,:);
        siz1 = [i1-i0+1,siz(2),siz(3)];
    elseif dimId==2
        A_bloc = A(:,i01:i11,:);
        siz1 = [siz(1),i1-i0+1,siz(3)];
    elseif dimId==3
        A_bloc = A(:,:,i01:i11);
        siz1 = [siz(1),siz(2),i1-i0+1];
    end
    
    % compute feature parameters in the image bloc
    [p_bloc] = featureParam(A_bloc,r,str_out);
    
    
    % allocate memory to feature parameters at the first bloc
    if ibloc==1
        param = cell(np,3);
        ndim = zeros(np,1);
        for i=1:np
            %dimension of the feature parameter
            ndim(i) = size(p_bloc{i},2);
            %allocate memory according to dimension
            for j=1:ndim(i)
                param{i,j} = zeros(siz,'single');
            end
        end
    end
    
    % place the results back to the corresponding positions
    id = false(size(A_bloc));
    if dimId==1
        id(m0:end-m1,:,:) = 1;
        for i=1:np
            for j=1:ndim(i)
                param{i,j}(i0:i1,:,:) = reshape(p_bloc{i}(id,j),siz1);
            end
        end
    elseif dimId==2
        id(:,m0:end-m1,:) = 1;
        for i=1:np
            for j=1:ndim(i)
                param{i,j}(:,i0:i1,:) = reshape(p_bloc{i}(id(:),j),siz1);
            end
        end
    elseif dimId==3
        id(:,:,m0:end-m1) = 1;
        for i=1:np
            for j=1:ndim(i)
                param{i,j}(:,:,i0:i1) = reshape(p_bloc{i}(id(:),j),siz1);
            end
        end
    end

    if over==1; break; end
end


% ouputs
if nargout==1
    varargout{1} = param;
else
    for i=1:np
        if ndim(i)>1
            varargout{i} = param(i,:);
        elseif ndim(i)==1
            varargout{i} = param{i,1};
        end
    end
end

