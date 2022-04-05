function [] = rmVoxReg3D(bw)
%remove connected regions by mouse in a 3D binary image
%
%to be continued...

% connected-component labeling
CC = bwconncomp(bw,26);
lm = labelmatrix(CC);

% choose polygon vertices
x = [];
y = [];
z = [];
figure;
for i=SLlst
    SL = A(:,:,i);
    imshow(SL);title(['slice number: ',num2str(i)]);
    
    w=waitforbuttonpress;
    if w==1; continue; end %keypress -> next slice
    
    [xi,yi] = ginput;
    zi = zeros(size(xi)) + i;
    x = [x;xi];
    y = [y;yi];
    z = [z;zi];
end
