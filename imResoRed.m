function [A] = imResoRed(A,w,str)
% This funciton will help in Corasegrain the image, with a special choice in gray value selection
%reduce 3D image resolution by 
%   1 - finding majority within a local neighborhood
%   2 - averaging withing a local neighborhood
%
%   [A] = imResoRed(A,w,str)
%   ------------------------
%
%   Inputs:
%       >  A : image to be processed
%       >  w : window size, vector of 1x3
%       > str: string specifying the method for resolution reduction
%              - 'majority': find the majority within each neighborhood
%              - 'average' : average within each neighborhood
%
%   Outputs:
%       > A : image of reduced resolution
%
% Keerthi Krishna PARVATHANENI 2018.04.11
%

siz = size(A);

% w(i) must be odd
for i=1:3
    if mod(w(i),2)==0;  w(i)=w(i)+1; end
end

% find the center location of each local neighborhood
nw = floor( siz ./ w ); %number of window in every dimension (new image resolution)

[cY,cX,cZ] = ndgrid( (w(1)+1)/2:w(1):nw(1)*w(1), ...
                     (w(2)+1)/2:w(2):nw(2)*w(2), ...
                     (w(3)+1)/2:w(3):nw(3)*w(3) );
cInd = sub2ind(siz,cY(:),cX(:),cZ(:));

% find the neighborhood indices
nbhoodInd = ind2indneighb(siz,cInd,'cube',(w-1)./2);


if strcmp(str,'majority')
    % majority within neighborhood
    A = mode(A(nbhoodInd),2);
    A = reshape(A,nw);
elseif strcmp(str,'average')
    % averaging within neighborhood
    A = mean(A(nbhoodInd),2);
    A = reshape(A,nw);
else
    error('string must be ''majority'' or ''average''')
end
