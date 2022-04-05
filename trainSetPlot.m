function trainSetPlot(Tsets,A)
%visualize the selected training set
%
%   trainSetPlot(Tsets,A)
%   ---------------------
%
%   Inputs:
%       > Tsets : coordinates of training sets (see trainSets.m)
%       >     A : Image in which the training sets are defined
%
%   Ouputs: figures
%
% Keerthi Krishna PARVATHANENI 2018.02.26
%

ncl = size(Tsets,1); %nomber of sets/clusters

for icl=1:ncl
    % plot the first slice
    isl = Tsets(icl,3);
    minA = min(min(A(:,:,isl)));
    maxA = max(max(A(:,:,isl)));
    figure;
    subplot(2,1,1);imshow(A(:,:,isl),[minA,maxA]);hold on;
    rectangle('Position',[Tsets(icl,1:2),(Tsets(icl,4:5)-Tsets(icl,1:2))],...
              'EdgeColor','r')
    % plot the last slice
    isl = Tsets(icl,6);
    minA = min(min(A(:,:,isl)));
    maxA = max(max(A(:,:,isl)));
    
    subplot(2,1,2);imshow(A(:,:,isl),[minA,maxA]);hold on;
    rectangle('Position',[Tsets(icl,1:2),(Tsets(icl,4:5)-Tsets(icl,1:2))],...
              'EdgeColor','r')
    xlabel(['training set - ',num2str(icl)]);
end





