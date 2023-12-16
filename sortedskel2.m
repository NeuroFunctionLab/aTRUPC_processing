function [coor_AP,coor_disp] = sortedskel2(BW_ske)

ske_end = bwmorph(BW_ske,'endpoints');
figure,imshow(ske_end,[]); title('ends of skeleton')

% find anterior point
[row_end,col_end] = find(ske_end==1);
coor_end = [row_end,col_end];
[~,index] = min(col_end(:));
coor_endA = coor_end(index,:);
[row,col] = find(BW_ske==1); coor = [row,col];
% find the shortest A->P trace

[disp , index] = min(sum(abs(coor - repmat(coor_endA,[size(coor,1),1])),2));
coor_new = coor;
coor_newend = coor_end(1,:);
coor_AP = [];
coor_disp = [];

for i =1:length(row)
    coor_new(index,:) = [];
    coor_AP = cat(1,coor_AP,coor_newend);
    coor_disp = cat(1,coor_disp,disp);
    dist = abs(coor_new - repmat(coor_newend(1,:),[size(coor_new,1),1]));
    [disp, index] = min(sqrt(dist(:,1).^2 + dist(:,2).^2));
    coor_newend = coor_new(index,:);
end