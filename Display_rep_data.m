clear;

%% aTRUPC choose
% data_path = 'D:\jsong\COVID_data\20211021\aTRUPC_mid_sagittal_real_imag_Auto_V13\results_real_imag.mat'; % COVID
data_path = 'D:\jsong\COVID_data\Control\NewControl\R01_20221111_ND\rawdata\aTRUPC_mid_sagittal_real_imag_Auto_V13\results_real_imag.mat'; % control
load(data_path);

Im_caxis = [0, 500];
figure,imagesc(abs(mag1{1}));colormap gray;axis image;axis off;caxis(Im_caxis);    
% set (gca,'Position',[0 0 1 1]);                                      
% print(gcf,'Phase reference.jpg', '-djpeg', '-r300');
figure,imagesc(abs(mag2{1}));colormap gray;axis image;axis off;caxis(Im_caxis);


Im_caxis = [0, 100];
figure,imagesc(abs(CD_complex{1}));colormap gray;axis image;axis off;caxis(Im_caxis);
hold on;
for iROI = 1:length(BW)
    plot(ROIx{iROI}, ROIy{iROI}, 'r', 'linewidth', 2);
end
hold off;

row_range = [29:186];
col_range = [21:245];
for i_eTE = 1:length(CD_complex)
    figure,imagesc(abs(CD_complex{i_eTE}(row_range, col_range)));colormap gray;axis image;axis off;caxis(Im_caxis);
end


for i_eTE = 1:length(CD_complex)
    figure,imagesc(abs(CD_complex{i_eTE}));colormap gray;axis image;axis off;caxis(Im_caxis);
end
%% ROI selection

fig_ROI3 = figure('name','select ROIs');
    addToolbarExplorationButtons(fig_ROI3); % this is important, otherwise using zoom-in during roipoly() is super annoying!
%     imshow(mask_vessel_final);
    imagesc(abs(CD_complex{1}));colormap gray;axis image;axis off;caxis(Im_caxis);
    axis image;colormap gray;
    set(gca,'fontsize',15);
    axis off;
    
    nROI = 3;
    BW = cell(3,1);
    ROIx = cell(3,1);
    ROIy = cell(3,1);
    hold on;
    for iROI = 1:nROI
        if iROI == 1
            title(sprintf('draw ROI#%d / %d', iROI, 1));
        elseif iROI == 2
            title(sprintf('draw ROI#%d / %d', iROI, 2));
        else 
            title(sprintf('draw ROI#%d / %d', iROI, 3));    
        end    
        [BW{iROI}, ROIx{iROI}, ROIy{iROI}] = roipoly();
        plot(ROIx{iROI},ROIy{iROI},'r','linewidth',2);
%         text(mean(ROIx{iROI}), mean(ROIy{iROI}), sprintf('%d', iROI), 'color', [0.929, 0.694, 0.125],...
%             'fontsize', 20, 'fontweight', 'bold');
    end
    hold off;
    title(sprintf('ROIs'));set(gca,'fontsize',15);
    zoom out;




%% Display ROI fitting
ROI2disp = [1 3 6];

for iROI = 1:length(ROI2disp)
    tx = linspace(0, max(eTE(:))+20, 1000);
    sig_fit1 = sig_ROI_mag_avg(1, ROI2disp(iROI))*exp(-(tx-eTE(1))/T2_roi_AllRep(ROI2disp(iROI)));
    dot_color = [56, 92, 163]/255;
    figure,plot(eTE, sig_ROI_mag_avg(:, ROI2disp(iROI)), 'o', 'linewidth', 1.2, 'markerfacecolor', dot_color, 'markeredgecolor', dot_color, 'markersize', 20);
    hold on;plot(tx, sig_fit1, '--', 'linewidth', 3, 'color', [0,0,0]);
    set(gca, 'fontsize', 25, 'fontweight', 'bold', 'linewidth', 3, 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
    box off;
    hold off;
    set(gca,'xcolor',[0,0,0]);
    set(gca,'ycolor',[0,0,0]);
    xticks([0, 40, 80]);
%     ylim([0, 50]);
end


%% OEF map
Ya = 98;
OEFmap = (Ya - Yvmap_corrected)./Ya*100;
OEFmap = OEFmap.*double(Yvmap_corrected > 0);
figure,imagesc(OEFmap);
OEFcaxis = [0,100];
OEFcmap = jet;
ind_0 = ceil(size(OEFcmap,1)*(0-OEFcaxis(1))/(OEFcaxis(2)-OEFcaxis(1))); % row index in the colormap for value 0
ind_0 = min(max(ind_0, 1),size(OEFcmap,1));
OEFcmap(ind_0,:) = [0,0,0]; % set color for value 0 to black
imagesc(OEFmap);axis image;
colorbar;set(gca,'fontsize',20,'fontweight','bold');
caxis(OEFcaxis);
colormap(OEFcmap);
axis off;

    

