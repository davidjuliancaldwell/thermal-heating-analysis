function [xx,yy,T_grid,x,y,T] = importfile_quickfield(filename)
% DJC - BLBT Summer Exchange 2017. This is a function to convert a 
% quickfield heating simulation file to one that can be analyzed in matlab
% e.g. Sequencefile01 = convert_demo('Sequencefile_01.csv');


%% Open the text file.
fileID = fopen(filename,'r');

startRow = 1;
header = cell2mat(textscan(fileID, '%f %f %f',1));

dataArray = textscan(fileID,'%f %f %f','HeaderLines',2,'ReturnOnError', false);

fclose(fileID);

x = dataArray{1};
y = dataArray{2};
T = dataArray{3};

dx = diff(x);
dx = (dx(dx>0));
dx = dx(1);
%total_x = max(x)-min(x);
total_x = header(1)*dx;
x_vec = [0:dx:total_x];
x_vec = x_vec(1:end-1);

dy = diff(y);
dy = (dy(dy>0));
dy = dy(1);
%total_y = max(y) - min(y);
total_y = header(2)*dx;
y_vec = [0:dy:total_y];
y_vec = y_vec(1:end-1);

[xx,yy] = meshgrid(x_vec,y_vec);

length_x = length(x_vec);
length_y = length(y_vec);

T = T -273;
T_grid = reshape(T,[length_x,length_y])';
figure;
imagesc(T_grid);
gcf;
set(gca,'YDir','reverse');
colorbar;
title('Simulated FEM Heating');


end