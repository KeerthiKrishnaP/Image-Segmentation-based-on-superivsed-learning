function [v0,COV] = covMatrix(v)
%Compute the centroid location and the covariance matrix of a set of
%feature parameters
%
%   [v0,COV] = covMatrix(v)
%   ------------------------
%   Inputs:
%       >   v: feature vector to be processed
%
%   Outputs:
%       >  v0: centroid locations
%       > COV: covariance matrix
%
%   Keerthi Krishna PARVATHANENI 2018.01.25
%

np = size(v,2); %number of feature parameters
nvx = size(v,1); %number of points (voxels)

% compute the centroid for each training cluster
v0 = mean(v,1);

% compute the covariance matrix for each training cluster
v1 = zeros(nvx,np,'single');
for ip = 1:np
    v1(:,ip) = v(:,ip)-v0(ip);
end
covmatrix = zeros(nvx,np,np,'single');
for ip = 1:np
    for jp = 1:np
        covmatrix(:,ip,jp) = v1(:,ip).*v1(:,jp);
    end
end
COV = sum(covmatrix,1) ./ (nvx-1); %why N-1,not just N?
COV = squeeze(COV);
