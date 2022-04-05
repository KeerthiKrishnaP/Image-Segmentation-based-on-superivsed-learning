function [idx_neighb] = ind2indneighb_p(varargin)
%This function returns the voxel list of the input points neighbors (only
%for 3D application, 2D not yet availabel)  
%
%   [idx_neighb] = ind2indneighb(siz,idx)
%   [idx_neighb] = ind2indneighb(siz,idx,shape,config)
%
%   Inputs:
%       siz: the size of image reference (ex. siz=size(VOL0))
%       idx: the indices of input points
%       shape: string defining the box shape {'cube'(defaut, with a
%              dimension of 11x11x11),'ellipsoid','cylinder'} 
%       config: vector defining the shape parameters ([row,col,band])
%               Attention: the dimensions for the cubic box are
%               2*config+1 !!!
%       periodicity: integer number (1,2,3) noting that the i-th direction
%                    is considered as periodic
%
%   Outputs:
%       idx_neighb: a matrix of mxn containing the voxel list of neighbors
%                   for all input points. m=length(idx), n=nb_vox_neighbor
%
% Keerthi Krishna PARVATHANENI  28/09/2015
% Keerthi Krishna PARVATHANENI  22/11/2016
%

[siz,idx,shape,config,periodicity] = ParseInputs(varargin{:});
% 
%disp(['neighborhood recognizing: ',shape,' ',num2str(2.*config+1)])

% type = ['uint',num2str(ceil(log2(prod(siz))))];
type = ceil(log2(prod(siz)));
if type<=8;
    type = 'uint8';
elseif type>8 && type<=16
    type = 'uint16';
elseif type>16 && type<=32
    type = 'uint32';
elseif type>32 && type<64
    type = 'uint64';
else
    disp 'the neighboring table is too heavy !!'
end
    
        [idxR,idxC,idxB] = ind2sub(siz,idx);

switch shape
    case 'cube'
        sizneighb = [config(1)*2+1,config(2)*2+1,config(3)*2+1];
        idx_neighb = zeros(length(idx),prod(sizneighb),type);
        sizcum = cumprod(siz);
        for i=1:prod(sizneighb)
            [ia,ib,ic] = ind2sub(sizneighb,i);
            ia = ia-config(1)-1;
            ib = ib-config(2)-1;
            ic = ic-config(3)-1;
            idx_neighb(:,i) = idx + ic*sizcum(2) + ib*sizcum(1) + ia;
            
            % out of limit verification
            outoflimitR = idxR+ia<=0 | idxR+ia>siz(1);
            outoflimitC = idxC+ib<=0 | idxC+ib>siz(2);
            outoflimitB = idxB+ic<=0 | idxB+ic>siz(3);
            outoflimit = outoflimitR | outoflimitC | outoflimitB;
            idx_neighb(outoflimit,i) = 0;
            
            % periodicity consideration
            if nnz(periodicity)==1 % only one periodicity
                
                if periodicity==1
                    i0 = idxC+ib>0 & idxC+ib<=siz(2) & ...
                         idxB+ic>0 & idxB+ic<=siz(3) ;
                    i00 = idxR+ia<=0 & i0;
                    idx_neighb(i00,i) = sub2ind( siz,...
                                                 idxR(i00)+ia+siz(1),...
                                                 idxC(i00)+ib,...
                                                 idxB(i00)+ic );
                    i00 = idxR+ia>siz(1) & i0;
                    idx_neighb(i00,i) = sub2ind( siz,...
                                                 idxR(i00)+ia-siz(1),...
                                                 idxC(i00)+ib,...
                                                 idxB(i00)+ic );
                elseif periodicity==2
                    i0 = idxR+ia>0 & idxR+ia<=siz(1) & ...
                         idxB+ic>0 & idxB+ic<=siz(3) ;
                    i00 = idxC+ib<=0 & i0;
                    idx_neighb(i00,i) = sub2ind( siz,...
                                                 idxR(i00)+ia,...
                                                 idxC(i00)+ib+siz(2),...
                                                 idxB(i00)+ic );
                    i00 = idxC+ib>siz(2) & i0;
                    idx_neighb(i00,i) = sub2ind( siz,...
                                                 idxR(i00)+ia,...
                                                 idxC(i00)+ib-siz(2),...
                                                 idxB(i00)+ic );
                elseif periodicity==3
                    i0 = idxR+ia>0 & idxR+ia<=siz(1) & ...
                         idxC+ib>0 & idxC+ib<=siz(2) ;
                    i00 = idxB+ic<=0 & i0;
                    idx_neighb(i00,i) = sub2ind( siz,...
                                                 idxR(i00)+ia,...
                                                 idxC(i00)+ib,...
                                                 idxB(i00)+ic+siz(3) );
                    i00 = idxB+ic>siz(3) & i0;
                    idx_neighb(i00,i) = sub2ind( siz,...
                                                 idxR(i00)+ia,...
                                                 idxC(i00)+ib,...
                                                 idxB(i00)+ic-siz(3) );
                end
                
            elseif nnz(periodicity)==2 % two periodicities
                
                if periodicity(1)==1 && periodicity(2)==2
                    idxR_neigh = idxR+ia;
                    idxC_neigh = idxC+ib;
                    idxB_neigh = idxB+ic;
                    
                    i0 = idxR_neigh<=0;
                    idxR_neigh(i0) = idxR_neigh(i0) + siz(1);
                    i0 = idxR_neigh>siz(1);
                    idxR_neigh(i0) = idxR_neigh(i0) - siz(1);
                    
                    i0 = idxC_neigh<=0;
                    idxC_neigh(i0) = idxC_neigh(i0) + siz(2);
                    i0 = idxC_neigh>siz(2);
                    idxC_neigh(i0) = idxC_neigh(i0) - siz(2);
                    
                    i0 = idxB+ic>0 & idxB+ic<=siz(3) ;

                    idx_neighb(i0,i) = sub2ind( siz,...
                                               idxR_neigh(i0),...
                                               idxC_neigh(i0),...
                                               idxB_neigh(i0) );
                    
                elseif periodicity(1)==1 && periodicity(2)==3
                    idxR_neigh = idxR+ia;
                    idxC_neigh = idxC+ib;
                    idxB_neigh = idxB+ic;
                    
                    i0 = idxR_neigh<=0;
                    idxR_neigh(i0) = idxR_neigh(i0) + siz(1);
                    i0 = idxR_neigh>siz(1);
                    idxR_neigh(i0) = idxR_neigh(i0) - siz(1);
                    
                    i0 = idxB_neigh<=0;
                    idxB_neigh(i0) = idxB_neigh(i0) + siz(3);
                    i0 = idxB_neigh>siz(3);
                    idxB_neigh(i0) = idxB_neigh(i0) - siz(3);
                    
                    i0 = idxC+ib>0 & idxC+ib<=siz(2) ;

                    idx_neighb(i0,i) = sub2ind( siz,...
                                               idxR_neigh(i0),...
                                               idxC_neigh(i0),...
                                               idxB_neigh(i0) );
                                           
                
                elseif periodicity(1)==2 && periodicity(2)==3
                    idxR_neigh = idxR+ia;
                    idxC_neigh = idxC+ib;
                    idxB_neigh = idxB+ic;
                    
                    i0 = idxC_neigh<=0;
                    idxC_neigh(i0) = idxC_neigh(i0) + siz(2);
                    i0 = idxC_neigh>siz(2);
                    idxC_neigh(i0) = idxC_neigh(i0) - siz(2);
                    
                    i0 = idxB_neigh<=0;
                    idxB_neigh(i0) = idxB_neigh(i0) + siz(3);
                    i0 = idxB_neigh>siz(3);
                    idxB_neigh(i0) = idxB_neigh(i0) - siz(3);
                    
                    i0 = idxR+ia>0 & idxR+ia<=siz(1) ;

                    idx_neighb(i0,i) = sub2ind( siz,...
                                               idxR_neigh(i0),...
                                               idxC_neigh(i0),...
                                               idxB_neigh(i0) );
                
                end
                
            elseif nnz(periodicity)==3 % three periodicities
                    idxR_neigh = idxR+ia;
                    idxC_neigh = idxC+ib;
                    idxB_neigh = idxB+ic;
                    
                    i0 = idxR_neigh<=0;
                    idxR_neigh(i0) = idxR_neigh(i0) + siz(1);
                    i0 = idxR_neigh>siz(1);
                    idxR_neigh(i0) = idxR_neigh(i0) - siz(1);
                    
                    i0 = idxC_neigh<=0;
                    idxC_neigh(i0) = idxC_neigh(i0) + siz(2);
                    i0 = idxC_neigh>siz(2);
                    idxC_neigh(i0) = idxC_neigh(i0) - siz(2);
                    
                    i0 = idxB_neigh<=0;
                    idxB_neigh(i0) = idxB_neigh(i0) + siz(3);
                    i0 = idxB_neigh>siz(3);
                    idxB_neigh(i0) = idxB_neigh(i0) - siz(3);
                    
                    idx_neighb(:,i) = sub2ind( siz,...
                                               idxR_neigh(:),...
                                               idxC_neigh(:),...
                                               idxB_neigh(:) );

            end
            
        end
    case 'ellipsoid'
        disp '"ELLIPSOID" shape box not yet available'
    case 'cylinder'
        disp '"CYLINDER" shape box not yet available'
end






%%%
%%% ParseInputs
%%%
function [siz,idx,shape,config,periodicity] = ParseInputs(varargin)

% Check the number of input arguments.
narginchk(2,5);

% Determine the shape and its configuration from the user supplied string
% and values. 
switch nargin
    case 2 %defaut
        siz = varargin{1};
        idx = varargin{2};
        shape = 'cube';
        config = [5 5 5];
        periodicity = [];
    case 4
        siz = varargin{1};
        idx = varargin{2};
        shape = varargin{3};
        config = varargin{4};
        periodicity = [];
    case 5
        siz = varargin{1};
        idx = varargin{2};
        shape = varargin{3};
        config = varargin{4};
        periodicity = varargin{5};
        periodicity = unique(periodicity);
end
shape = validatestring(shape,{'cube','ellipsoid','cylinder'},mfilename,'SHAPE',1);

