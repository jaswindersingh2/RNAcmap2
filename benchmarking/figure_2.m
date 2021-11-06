clc; clear all; close all;


addpath(strcat(pwd, '/DrosteEffect-CubeHelix-80f7824'));
addpath(strcat(pwd, '/Violinplot-Matlab-0.1'));

map = cubehelix(18);
map = map(1:13,:);
rng(2)
map = map(randperm(13),:);
map(1,:) =  [1 0 0]; 
map(3,:) = [0,0.5,0]; 


[ha, pos] = tight_subplot(1,1,[.1 .05],[.08 .12],[.05 .01]);
axes(ha(1));


data = readtable(strcat('figure_2_data/blastn_no_hit.csv'));
blastn = table2array(data(:,'F1_all_bps')); 

data = readtable(strcat('figure_2_data/direct_infernal_no_hit.csv'));
infernal = table2array(data(:,'F1_all_bps')); 

data = readtable(strcat('figure_2_data/RNAcmap_no_hit.csv'));
RNAcmap = table2array(data(:,'F1_all_bps')); 

data = readtable(strcat('figure_2_data/RNAcmap_meta_no_hit.csv'));
RNAcmap_meta = table2array(data(:,'F1_all_bps')); 

data = readtable(strcat('figure_2_data/RNAcmap2_meta_no_hit.csv'));
RNAcmap2_meta = table2array(data(:,'F1_all_bps')); 

h1 = Violin(blastn,0.5,'Width',0.2,'ViolinColor',map(5,:),'ShowMean',true);hold on;
h2 = Violin(infernal,1,'Width',0.2,'ViolinColor',map(4,:),'ShowMean',true);hold on;
h3 = Violin(RNAcmap,1.5,'Width',0.2,'ViolinColor',map(3,:),'ShowMean',true);
h4 = Violin(RNAcmap_meta,2.0,'Width',0.2,'ViolinColor',map(2,:),'ShowMean',true);
h5 = Violin(RNAcmap2_meta,2.5,'Width',0.2,'ViolinColor',map(1,:),'ShowMean',true);


data = readtable(strcat('figure_2_data/blastn_low.csv'));
blastn = table2array(data(:,'F1_all_bps')); 

data = readtable(strcat('figure_2_data/direct_infernal_low.csv'));
infernal = table2array(data(:,'F1_all_bps')); 

data = readtable(strcat('figure_2_data/RNAcmap_low.csv'));
RNAcmap = table2array(data(:,'F1_all_bps')); 

data = readtable(strcat('figure_2_data/RNAcmap_meta_low.csv'));
RNAcmap_meta = table2array(data(:,'F1_all_bps')); 

data = readtable(strcat('figure_2_data/RNAcmap2_meta_low.csv'));
RNAcmap2_meta = table2array(data(:,'F1_all_bps')); 

hold on;
h1 = Violin(blastn,3.5,'Width',0.2,'ViolinColor',map(5,:),'ShowMean',true);hold on;
h2 = Violin(infernal,4,'Width',0.2,'ViolinColor',map(4,:),'ShowMean',true);hold on;
h3 = Violin(RNAcmap,4.5,'Width',0.2,'ViolinColor',map(3,:),'ShowMean',true);
h4 = Violin(RNAcmap_meta,5.0,'Width',0.2,'ViolinColor',map(2,:),'ShowMean',true);
h5 = Violin(RNAcmap2_meta,5.5,'Width',0.2,'ViolinColor',map(1,:),'ShowMean',true);


data = readtable(strcat('figure_2_data/blastn_medium.csv'));
blastn = table2array(data(:,'F1_all_bps')); 

data = readtable(strcat('figure_2_data/direct_infernal_medium.csv'));
infernal = table2array(data(:,'F1_all_bps')); 

data = readtable(strcat('figure_2_data/RNAcmap_medium.csv'));
RNAcmap = table2array(data(:,'F1_all_bps')); 

data = readtable(strcat('figure_2_data/RNAcmap_meta_medium.csv'));
RNAcmap_meta = table2array(data(:,'F1_all_bps')); 

data = readtable(strcat('figure_2_data/RNAcmap2_meta_medium.csv'));
RNAcmap2_meta = table2array(data(:,'F1_all_bps')); 

hold on;
h1 = Violin(blastn,6.5,'Width',0.2,'ViolinColor',map(5,:),'ShowMean',true);hold on;
h2 = Violin(infernal,7,'Width',0.2,'ViolinColor',map(4,:),'ShowMean',true);hold on;
h3 = Violin(RNAcmap,7.5,'Width',0.2,'ViolinColor',map(3,:),'ShowMean',true);
h4 = Violin(RNAcmap_meta,8.0,'Width',0.2,'ViolinColor',map(2,:),'ShowMean',true);
h5 = Violin(RNAcmap2_meta,8.5,'Width',0.2,'ViolinColor',map(1,:),'ShowMean',true);



data = readtable(strcat('figure_2_data/blastn_high.csv'));
blastn = table2array(data(:,'F1_all_bps')); 

data = readtable(strcat('figure_2_data/direct_infernal_high.csv'));
infernal = table2array(data(:,'F1_all_bps')); 

data = readtable(strcat('figure_2_data/RNAcmap_high.csv'));
RNAcmap = table2array(data(:,'F1_all_bps')); 

data = readtable(strcat('figure_2_data/RNAcmap_meta_high.csv'));
RNAcmap_meta = table2array(data(:,'F1_all_bps')); 

data = readtable(strcat('figure_2_data/RNAcmap2_meta_high.csv'));
RNAcmap2_meta = table2array(data(:,'F1_all_bps')); 

hold on;
h1 = Violin(blastn,9.5,'Width',0.2,'ViolinColor',map(5,:),'ShowMean',true);hold on;
h2 = Violin(infernal,10,'Width',0.2,'ViolinColor',map(4,:),'ShowMean',true);hold on;
h3 = Violin(RNAcmap,10.5,'Width',0.2,'ViolinColor',map(3,:),'ShowMean',true);
h4 = Violin(RNAcmap_meta,11.0,'Width',0.2,'ViolinColor',map(2,:),'ShowMean',true);
h5 = Violin(RNAcmap2_meta,11.5,'Width',0.2,'ViolinColor',map(1,:),'ShowMean',true);

ax = gca;
ax.FontSize = 12; 

xticks([0 1.5 4.5 7.5 10.5]);
xticklabels({'0','No hit RNAs','Low Neff RNAs', 'Medium Neff RNAs', 'High Neff RNAs','10'});
% xlabel("Test sets")

yticks([0 0.2 0.4 0.6 0.8 1.0]);
yticklabels({'0','0.2','0.4','0.6','0.8','1.0'});
ylabel("F1-score");

grid on;

annotation('textbox', [0.20 0.95 0.04 0.025], 'String',{''}, 'FontSize',12, 'FontName','Arial', 'EdgeColor',[1 1 1], 'LineWidth',1,...
    'BackgroundColor',map(5,:), 'Color',[1 0 1], 'FaceAlpha', 0.6);
annotation('textbox', [0.35 0.95 0.04 0.025], 'String',{''}, 'FontSize',12, 'FontName','Arial', 'EdgeColor',[1 1 1], 'LineWidth',1,...
    'BackgroundColor',map(4,:), 'Color', [0.3635    0.4788    0.1841], 'FaceAlpha', 0.6);
annotation('textbox', [0.50 0.95 0.04 0.025], 'String',{''}, 'FontSize',12, 'FontName','Arial', 'EdgeColor',[1 1 1], 'LineWidth',1,...
    'BackgroundColor',map(3,:), 'Color', [0.3635    0.4788    0.1841], 'FaceAlpha', 0.6);
annotation('textbox', [0.65 0.95 0.04 0.025], 'String',{''}, 'FontSize',12, 'FontName','Arial', 'EdgeColor',[1 1 1], 'LineWidth',1,...
    'BackgroundColor',map(2,:), 'Color', [0.3635    0.4788    0.1841], 'FaceAlpha', 0.6);
annotation('textbox', [0.80 0.95 0.04 0.025], 'String',{''}, 'FontSize',12, 'FontName','Arial', 'EdgeColor',[1 1 1], 'LineWidth',1,...
    'BackgroundColor',map(1,:), 'Color', [0.3635    0.4788    0.1841], 'FaceAlpha', 0.6);

annotation('textbox', [0.19 0.92 0.4 0.02], 'String',{'BLAST-N'}, 'FontSize',12, 'FontName','Arial', 'EdgeColor',[1 1 1], 'LineWidth',1,...
    'BackgroundColor',[1  1 1], 'Color',[0 0 0]);
annotation('textbox', [0.34 0.92 0.04 0.02], 'String',{'INFERNAL'}, 'FontSize',12, 'FontName','Arial', 'EdgeColor',[1 1 1], 'LineWidth',1,...
    'BackgroundColor',[1  1 1], 'Color', [0 0 0]);
annotation('textbox', [0.49 0.92 0.04 0.02], 'String',{'RNAcmap'}, 'FontSize',12, 'FontName','Arial', 'EdgeColor',[1 1 1], 'LineWidth',1,...
    'BackgroundColor',[1  1 1], 'Color', [0 0 0]);
annotation('textbox', [0.64 0.92 0.04 0.02], 'String',{'RNAcmap*'}, 'FontSize',12, 'FontName','Arial', 'EdgeColor',[1 1 1], 'LineWidth',1,...
    'BackgroundColor',[1  1 1], 'Color', [0 0 0]);
annotation('textbox', [0.79 0.92 0.04 0.02], 'String',{'RNAcmap2'}, 'FontSize',12, 'FontName','Arial', 'EdgeColor',[1 1 1], 'LineWidth',1,...
    'BackgroundColor',[1  1 1], 'Color', [0 0 0]);

set(gcf, 'Position', [800, 600, 1500, 500]);
% saveas(gcf, 'figure_2_y', 'epsc');
saveas(gcf, 'figure_2', 'png'); % .png save because facealpha in annotation don't save in .eps format.

