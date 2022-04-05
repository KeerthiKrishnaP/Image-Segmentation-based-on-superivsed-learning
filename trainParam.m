function [v0T,COVT] = trainParam(T,r,str)
%Train the centroid and covariance matrix using a set of training clusters
%
%   [v0T,COVT] = trainParam(T,r,str_vec)
%   ------------------------------------
%   Inputs:
%       >    T: cell array, containing the training clusters (images)
%       >    r: scalar, specifying the window "radius" for structure tensor 
%       >  str: string array, specifying the parameters in the feature
%               vector 
%
%   Output:
%       >  v0T: centoids of the training clusters
%       > COVT: covariance matrice of the training clusters
%
%   Keerthi Krishna PARVATHANENI 2018.0126
%

disp(['--> Training with a vector: ',strjoin(str)])
ncl = length(T); %number of traning clusters
np = length(str); %number of parameters in the feature vector

v0T = zeros(ncl,np);
COVT = zeros(ncl,np,np);
for icl=1:ncl
    disp([blanks(4),num2str(icl),'-th training cluster:'])
    % build the feature vector for each training cluster
    sizT = size(T{icl});
    sizW = 2.*[r r r]+1;
    id = sizT-sizW<0; %find the direction where the size is smaller than the window side
    if id(1)>0
        T{icl} = cat(1,flip(T{icl}(1:r,:,:),1),T{icl},...
                       flip(T{icl}(end-r:end,:,:),1));
    end
    if id(2)>0
        T{icl} = cat(2,flip(T{icl}(:,1:r,:),2),T{icl},...
                       flip(T{icl}(:,end-r:end,:),2));
    end
    if id(3)>0
        T{icl} = cat(3,flip(T{icl}(:,:,1:r),3),T{icl},...
                       flip(T{icl}(:,:,end-r:end),3));
    end
    sizT = size(T{icl});
    
    vT = featureParam(T{icl},r,str);
    
    for ip=1:np %remove the image borders
        vT{ip} = reshape(vT{ip},sizT);
        vT{ip}=vT{ip}(r+2:end-r,r+2:end-r,r+2:end-r);
        vT{ip} = vT{ip}(:);
    end
    
    vT = cell2mat(vT');
    
    % centroid location & covariance matrix
    [v0_icl,COV_icl] = covMatrix(vT);
    v0T(icl,:) = v0_icl;
    COVT(icl,:,:) = COV_icl;
    
    % display
    disp([blanks(4),'centroid: ',num2str(v0T(icl,:))]);
    disp([blanks(4),'covmatrix: ']);
    disp(squeeze(COVT(icl,:,:)))
end
