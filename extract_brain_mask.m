function BW_brain = extract_brain_mask(varargin)
% Jie Song
% 09/20/2023

% bextract=1 manually selection
% bextract=0 automatic extraction, mophological operation
Image = varargin{1};

if nargin > 1
    bextract = varargin{2};
else
    bextract = 0;
end

if bextract
    fig_brain_ROI = figure('name', 'Skull Stripping');
    imshow(Image, []);hold on;
    addToolbarExplorationButtons(fig_brain_ROI); % this is important, otherwise using zoom-in during roipoly() is super annoying!

    NumROI = 2;                       
    BW_brain_drawn = cell(NumROI, 1);
    x_brain = cell(NumROI, 1);
    y_brain = cell(NumROI, 1);
    title(sprintf('Draw brain mask: #%d', 1));set(gca,'fontsize',15);
    [BW_brain_drawn{1}, x_brain{1}, y_brain{1}] = roipoly();
    plot(x_brain{1}, y_brain{1}, 'r', 'linewidth', 2); 
    hold off;
    drawnow;
    zoom out;
    
% Thresholding
    brain_thresh = median(Image(BW_brain_drawn{1}))*0.5;
    b = Image.*BW_brain_drawn{1} > brain_thresh;
    BW_brain=imclose(b,strel('disk',18)); 
    BW_brain = BW_brain_drawn{1};
else
    % segment the brain mask 
    Im_norm = Image./max(Image(:)); 
    k1=Im_norm > graythresh(Im_norm); 
    SE=strel('disk',7,4);
    k2=imopen(k1,SE); 
    [b,NUM]=bwlabel(k2);% apply connected component analysis.
    for i_b = 1:NUM
        b_num(i_b) = length(find(b(:)==i_b));
    end
    [~,index] = max(b_num);
    b(b~=index)=0;
    b(b==index)=1;
    BW_brain=imclose(b,strel('disk',18)); 
end 