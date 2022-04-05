function [varargout] = featureParam(A,r,str)
%Compute the feature parameters of a 3-D image
%
%   [varargout] = featureParam(A,r,str)
%   -----------------------------------
%   Inputs:
%       >   A: 3-D image to be processed
%       >   r: scalar, window "radius", the cubic windows has a size of
%              (2*r+1)^3
%       > str: string cell, specifying the outputs. The order of the names
%              in str must be the same as those of outputs.
%              You can choose among:
%              'S' - structure tensor
%              'AVG' - average gray level
%              'beta' - anisotropy parametery
%              'lambda' - eigenvalues with lambda1<lambda2<lambda3
%              'vec1' - eigenvector corresponding to lambda1
%              'vec2' - eigenvector corresponding to lambda2
%              'vec3' - eigenvector corresponding to lambda3
%              'phi' - orientation angle about x-axis (projected in XY)
%              'psi' - orientation angle about z-axis (3D)
%
%   Outputs:
%       see Inputs -> str
%       note: if nargout==1, but length(str)>1, the parameters will be
%             saved into a cell array.
%
%   Keerthi Krishna PARVATHANENI 2018.01.25
%

disp('----------------------------------');
% compute the structure tensor --------------------------------------------
disp('    structure tensor computing');
tic
S = structTensor(A,r);
toc

% compute the eigenthings and the anisotropy parameter --------------------
disp('    eigenthings computing');
tic
[lambda,vec1,vec2,vec3,beta] = eigTensor(S);
toc

% compute the average gray level ------------------------------------------
disp('    average gray level computing');
tic
avg = avgI(A,S.r);
avg = avg(:);
toc

% orientation angles ------------------------------------------------------
disp('    orientation angles computing');
tic
phi = atan(abs(vec1(:,2)./vec1(:,1)));
psi = acos(abs(vec1(:,3)));
toc

% outputs -----------------------------------------------------------------
iS = find(strcmpi(str,'S'));
iLAMBDA = find(strcmpi(str,'lambda'));
iVEC1 = find(strcmpi(str,'vec1'));
iVEC2 = find(strcmpi(str,'vec2'));
iVEC3 = find(strcmpi(str,'vec3'));
iBETA = find(strcmpi(str,'beta'));
iAVG = find(strcmpi(str,'AVG'));
iPHI = find(strcmpi(str,'phi'));
iPSI = find(strcmpi(str,'psi'));

if nargout==1
    outputs = cell(length(str),1);
    if ~isempty(iS);      outputs{iS} = S;           end
    if ~isempty(iLAMBDA); outputs{iLAMBDA} = lambda; end
    if ~isempty(iVEC1);   outputs{iVEC1} = vec1;     end
    if ~isempty(iVEC2);   outputs{iVEC2} = vec2;     end
    if ~isempty(iVEC3);   outputs{iVEC3} = vec3;     end
    if ~isempty(iBETA);   outputs{iBETA} = beta;     end
    if ~isempty(iAVG);    outputs{iAVG} = avg;       end
    if ~isempty(iPHI);    outputs{iPHI} = phi;       end
    if ~isempty(iPSI);    outputs{iPSI} = psi;       end
    varargout{1} = outputs;
elseif nargout>1
    if ~isempty(iS);      varargout{iS} = S;           end
    if ~isempty(iLAMBDA); varargout{iLAMBDA} = lambda; end
    if ~isempty(iVEC1);   varargout{iVEC1} = vec1;     end
    if ~isempty(iVEC2);   varargout{iVEC2} = vec2;     end
    if ~isempty(iVEC3);   varargout{iVEC3} = vec3;     end
    if ~isempty(iBETA);   varargout{iBETA} = beta;     end
    if ~isempty(iAVG);    varargout{iAVG} = avg;       end
    if ~isempty(iPHI);    varargout{iPHI} = phi;       end
    if ~isempty(iPSI);    varargout{iPSI} = psi;       end
end

