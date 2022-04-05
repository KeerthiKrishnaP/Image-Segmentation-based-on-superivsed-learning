function [P,covprod] = probFun(v,v0T,COVT)
%Compute the probability of a feature vector with respect to a reference
%cluster that is described by its centroid and covariance matrix
%
%   [P] = probFun(v,v0T,COVT)
%   -------------------------
%   Inputs:
%       >    v: feature vector
%       >  v0T: centroid of the reference cluster
%       > COVT: covariance matrix of the reference cluster
%
%   Outputs:
%       >    P: probability parameter
%
%   note: the formula is taken from [Straumit,2015], without the term of pi
%
%   Keerthi Krishna PARVATHANENI  2018.01.25
%

np = size(v,2); %number of parameters
nvx = size(v,1); %number of points (voxels)

% inverse the reference covariance matrix
COVT_inv = inv(COVT);

% probability function: ref. [Straumit,2015]
v1 = zeros(nvx,np,'single');
for ip = 1:np
    v1(:,ip) = v(:,ip)-v0T(ip);
end

covprod = zeros(nvx,1,'single');
for ip = 1:np
    for jp = 1:np
        covprod = covprod + v1(:,ip).*COVT_inv(ip,jp).*v1(:,jp);
    end
end

P = exp(covprod./(-2)) ./ sqrt( det(COVT) );
