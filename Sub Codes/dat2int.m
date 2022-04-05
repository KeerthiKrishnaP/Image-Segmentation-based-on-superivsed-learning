function dat = dat2int(varargin)
%set the data to an apropriate integer type to reduce the memory use
%
%   dat = dat2int(dat)
%   dat = dat2int(dat,minmax)
%   -----------------------
%
%   Inputs:
%       >    dat : input data
%       > minmax : (optional) the min/max values in the input data
%
%   Outputs:
%       >    dat : output data with the new integer type
%
% Keerthi Krishna PARVATHANENI  2018.02.23
%
narginchk(1,2);

dat = varargin{1};

if nargin==1
    minmax = [min(dat(:)), max(dat(:))];
elseif nargin==2
    minmax = varargin{2};
end


if minmax(1)<0 %int*
    if minmax(1)>=intmin('int8') && minmax(2)<=intmax('int8')
        dat = int8(dat);
    elseif minmax(1)>=intmin('int16') && minmax(2)<=intmax('int16')
        dat = int16(dat);
    elseif minmax(1)>=intmin('int32') && minmax(2)<=intmax('int32')
        dat = int32(dat);
    elseif minmax(1)>=intmin('int64') && minmax(2)<=intmax('int64')
        dat = int64(dat);
    end
elseif minmax(1)>=0 %uint*
    if minmax(2)<=intmax('uint8')
        dat = uint8(dat);
    elseif minmax(2)<=intmax('uint16')
        dat = uint16(dat);
    elseif minmax(2)<=intmax('uint32')
        dat = uint32(dat);
    elseif minmax(2)<=intmax('uint64')
        dat = uint64(dat);
    end
end

return
end