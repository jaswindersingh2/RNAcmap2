clc; clear all; close all;

addpath(strcat(pwd, '/DrosteEffect-CubeHelix-80f7824'));

map = cubehelix(40);
rng(2)
for k = 1:length(map)
    map(2*k-1,:) = [0,0.5,0]; 
    map(2*k,:) =  [1 0 0]; 
end

[ha, pos] = tight_subplot(3,1,[.1 .05],[.06 .1],[.03 .01]);


save_data = zeros(10, 2, 245);

for n = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    data1 = readtable(strcat('figure_3_data/RNAcmap_mfdca_top_L_by_', char(int2str(n)), '.csv'));
    data2 = readtable(strcat('figure_3_data/RNAcmap2_meta_mfdca_top_L_by_', char(int2str(n)), '.csv'));

    temp = [table2array(data1(:,'F1_all_bps')), table2array(data2(:,'F1_all_bps'))];
    save_data(n, :, :) = temp';

end


axes(ha(1));

y = [save_data];
x = 1:10;
b = boxplot2(y, x,'Notch','on');
grid on;
axis tight; 

text(-0.03, 1.2, char('\bfA'),'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize',14);

xticks([1 2 3 4 5 6 7 8 9 10]);
set(gca,'XTickLabel',{'L/1', 'L/2', 'L/3', 'L/4', 'L/5', 'L/6', 'L/7', 'L/8', 'L/9', 'L/10'});
xlabel('Top L/n contacts','FontSize',12);

ylabel('F1-Score','FontSize',12);
ylim([0 1]);
set(gca,'YTickLabel',{'0', '0.2', '0.4', '0.6', '0.8', '1.0'});

h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(end+1-j),'XData'),get(h(end+1-j),'YData'),map(j,:));
end


hold on
for j=1:length(h)
    xdat = get(h(end+1-j),'XData');
    ydat = get(h(end+1-j),'YData');
    handle=line([xdat(1),xdat(6)],[ydat(1),ydat(6)]);
    set(handle,'Color','k');
end


save_data = zeros(10, 2, 245);

for n = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    data1 = readtable(strcat('figure_3_data/RNAcmap_mfdca_top_L_by_', char(int2str(n)), '.csv'));
    data2 = readtable(strcat('figure_3_data/RNAcmap2_meta_mfdca_top_L_by_', char(int2str(n)), '.csv'));

    temp = [table2array(data1(:,'PR_all_bps')), table2array(data2(:,'PR_all_bps'))];
    save_data(n, :, :) = temp';

end

axes(ha(2));

y = save_data;
x = 1:10;
b = boxplot2(y, x,'Notch','on');
grid on;
axis tight; 

text(-0.03, 1.2, char('\bfB'),'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize',14);

xticks([0 1 2 3 4 5 6 7 8 9 10]);
set(gca,'XTickLabel',{' ', 'L/1', 'L/2', 'L/3', 'L/4', 'L/5', 'L/6', 'L/7', 'L/8', 'L/9', 'L/10'});
xlabel('Top L/n contacts','FontSize',12);

ylabel('Precision','FontSize',12);
ylim([0 1]);
set(gca,'YTickLabel',{'0', '0.2', '0.4', '0.6', '0.8', '1.0'});

h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(end+1-j),'XData'),get(h(end+1-j),'YData'),map(j,:));
end


hold on
for j=1:length(h)
    xdat = get(h(end+1-j),'XData');
    ydat = get(h(end+1-j),'YData');
    handle=line([xdat(1),xdat(6)],[ydat(1),ydat(6)]);
    set(handle,'Color','k');
end


save_data = zeros(10, 2, 245);

for n = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    data1 = readtable(strcat('figure_3_data/RNAcmap_mfdca_top_L_by_', char(int2str(n)), '.csv'));
    data2 = readtable(strcat('figure_3_data/RNAcmap2_meta_mfdca_top_L_by_', char(int2str(n)), '.csv'));

    temp = [table2array(data1(:,'SN_all_bps')), table2array(data2(:,'SN_all_bps'))];
    save_data(n, :, :) = temp';

end

axes(ha(3));

y = save_data;
x = 1:10;
b = boxplot2(y, x,'Notch','on');
grid on;
axis tight; 

text(-0.03, 1.2, char('\bfC'),'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize',14);

xticks([0 1 2 3 4 5 6 7 8 9 10]);
set(gca,'XTickLabel',{' ', 'L/1', 'L/2', 'L/3', 'L/4', 'L/5', 'L/6', 'L/7', 'L/8', 'L/9', 'L/10'});
xlabel('Top L/n contacts','FontSize',12);

ylabel('Sensitivity','FontSize',12);
ylim([0 1]);
set(gca,'YTickLabel',{'0', '0.2', '0.4', '0.6', '0.8', '1.0'});

h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(end+1-j),'XData'),get(h(end+1-j),'YData'),map(j,:));
end


hold on
for j=1:length(h)
    xdat = get(h(end+1-j),'XData');
    ydat = get(h(end+1-j),'YData');
    handle=line([xdat(1),xdat(6)],[ydat(1),ydat(6)]);
    set(handle,'Color','k');
end


annotation('textbox', [0.45 0.95 0.04 0.025], 'String',{''}, 'FontSize',12, 'FontName','Arial', 'EdgeColor',[1 1 1], 'LineWidth',1,...
    'BackgroundColor',[0,0.5,0], 'Color',[1 0 1]);
annotation('textbox', [0.55 0.95 0.04 0.025], 'String',{''}, 'FontSize',12, 'FontName','Arial', 'EdgeColor',[1 1 1], 'LineWidth',1,...
    'BackgroundColor',[1 0 0], 'Color', [0.3635    0.4788    0.1841]);

annotation('textbox', [0.44 0.92 0.4 0.02], 'String',{'RNAcmap'}, 'FontSize',12, 'FontName','Arial', 'EdgeColor',[1 1 1], 'LineWidth',1,...
    'BackgroundColor',[1  1 1], 'Color',[0 0 0]);
annotation('textbox', [0.54 0.92 0.04 0.02], 'String',{'RNAcmap2'}, 'FontSize',12, 'FontName','Arial', 'EdgeColor',[1 1 1], 'LineWidth',1,...
    'BackgroundColor',[1  1 1], 'Color', [0 0 0]);

set(gcf, 'Position', [800, 600, 1600, 800]);
% saveas(gcf, 'figure_3', 'png');
