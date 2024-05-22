function [out, out2, coor_AP,coor_disp]= sortedskel(varargin)
%(BW_SSS,preprocessind,branchlimit)
% Based on inbuilt function bwskel & bwmorph in imaging processing box
% sorted all the skeleton points A to P (left to right), tested on 
% SSS, Straight-sinus
% input: 
% BW_SSS: vessel mask or vessel skeleton 
% preprocessind: 1 for extracting the skeleton
% branchlimit: remove unwanted branches on main trace

% output:
% out: branched skeleton: vessel skeleton
% out2: non-branched skeleton
% coor_AP: sorted A->P points array of non-branched skeleton

% Jie Song
% 07/04/2023

BW_SSS = varargin{1};
preprocessind = varargin{2};
if nargin > 2
    branchlimit = varargin{3};
end

out = BW_SSS;
out2 = BW_SSS;
% Part 1 - pre-processing for skeleton extraction
if preprocessind == 1
    BW_SSS = bwareaopen(BW_SSS, 20);
    for i = 1:3
        BW_SSS = bwmorph(BW_SSS,'spur');
    end
    BW_SSS = logical(imfill(BW_SSS,'holes'));
    % extract skeleton
    out = bwskel(BW_SSS);
    out2 = bwskel(BW_SSS,'MinBranchLength',branchlimit); % the size of the disk needs modification
    % figure,imshow([out,out2],[]); title('skeleton of SSS')
end

out2_end = bwmorph(out2,'endpoints');
% figure,imshow([out2,out2_end],[]); title('ends of skeleton')

% Comparison
% out = bwskel(logical(BW_SSS));
% [tryx1,tryy1] = find(out);
% figure, imagesc(BW_SSS);axis off;axis image;colormap gray;hold on;
% plot(tryy1,tryx1,'.r','linewidth',2);
% hold off;
% 
% [tryx1,tryy1] = find(out2);
% figure, imagesc(BW_SSS);axis off;axis image;colormap gray;hold on;
% plot(tryy1,tryx1,'.r','linewidth',2);
% hold off;

% Part 2 - sort all the skeleton points: A to P
% find anterior point
[row_end,col_end] = find(out2_end==1);
coor_end = [row_end,col_end];
[~,index] = min(col_end(:));
coor_endA = coor_end(index,:);
[row,col] = find(out2==1); coor = [row,col];

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