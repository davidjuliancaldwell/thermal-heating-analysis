% 7-14-2017 - script to analyze heating profile in the thermal camera
% imaging data

%% load file
close all;clear all;clc

file = 9;

switch file
    case 1
        load('elec1.5d3.05p_fu_8.0mA_0001.mat'); % this was the one I started with
        yVals = [250:370];
        xVals = [130:295];
    case 2
        
        load('elec2.0d3.05p_fu_8.0mA_0001.mat');
        yVals = [240:380];
        xVals = [140:300];
    case 3
        load('grid2_optistim5_18return_26in_10in_1.75mm_0001.mat')
        yVals = [240:380];
        xVals = [140:270];
    case 4
        load('grid2_optistim3_18return_10_20_26in_1.75mm_0001.mat')
        yVals = [180:350];
        xVals = [140:300];
    case 5
        load('optistim4_dual_1.75_3.05_11_19_27_trial1_smallSection.mat')
        yVals = [1:size(stackedData,1)];
        xVals = [1:size(stackedData,2)-5];
    case 6
        load('C:\david\david\11-8-2017\2nd_gel_flipped_betweenGrids\2ndGel_betweenGrids_5_8_run1_converted\2ndGel_betweenGrids_5_8_run1.mat')
        xVals = [240:380];
        yVals = [180:300];
        % convert to centigrade
        stackedData = stackedData - 273.15;
    case 7
        load('C:\david\david\11-8-2017\2nd_gel_flipped_betweenGrids\2ndGel_betweenGrids_5_8_run2_converted\2ndGel_betweenGrids_5_8_run2.mat')
        xVals = [240:380];
        yVals = [180:300];
        % convert to centigrade
        stackedData = stackedData - 273.15;
    case 8
        load('C:\david\david\11-8-2017\2nd_gel_flipped_betweenGrids\2ndGel_betweenGrids_5_8_run3_converted\2ndGel_betweenGrids_5_8_run3.mat')
        xVals = [240:380];
        yVals = [180:300];
        % convert to centigrade
        stackedData = stackedData - 273.15;
    case 9
        xVals = [240:380];
        yVals = [180:300];
        load('C:\david\david\11-8-2017\2nd_gel_flipped_betweenGrids\2ndGel_betweenGrids_5_8_run1_converted\2ndGel_betweenGrids_5_8_run1.mat')
        stackedData = stackedData - 273.15;
        stackedDataSub1 = stackedData(xVals,yVals,:);
        clear stackedData;
        load('C:\david\david\11-8-2017\2nd_gel_flipped_betweenGrids\2ndGel_betweenGrids_5_8_run2_converted\2ndGel_betweenGrids_5_8_run2.mat')
        stackedData = stackedData - 273.15;
        stackedDataSub2 = stackedData(xVals,yVals,:);
        clear stackedData;
        load('C:\david\david\11-8-2017\2nd_gel_flipped_betweenGrids\2ndGel_betweenGrids_5_8_run3_converted\2ndGel_betweenGrids_5_8_run3.mat')
        stackedData = stackedData - 273.15;
        stackedDataSub3 = stackedData(xVals,yVals,:);
        clear stackedData;
        
end

stackedDataSub = stackedData(xVals,yVals,:);


smoothIt = 0;

% 2d smoothing
if smoothIt
    %h = 1/9*[0 1 0; 1 5 1; 0 1 0];
    h = 0.125*ones(3);
    for i = 1:size(stackedDataSub,3)
        stackedDataSub(:,:,i) = conv2(stackedDataSub(:,:,i),h,'same');
    end
    
    stackedDataSub = stackedDataSub(2:end-1,2:end-1,:);
end

plotIt = 1;

%% begin the analysis

% visualize
if plotIt
    
    fig1 = figure;
    for ind = 1:size(stackedData,3)
        
        subplot(3,1,1)
        imagesc(stackedData(:,:,ind));
        colorbar()
        drawnow()
        
        subplot(3,1,2)
        contour(stackedData(:,:,ind));
        set(gca,'YDir','reverse')
        colorbar()
        drawnow()
        
        subplot(3,1,3)
        histogram(reshape(stackedData(:,:,ind),[],1),100,'Normalization','probability')
        drawnow()
        
    end
end
%% look at subportion of thermal heating



if plotIt
    
    fig2 = figure;
    
    for ind = 1:size(stackedDataSub,3)
        
        subplot(3,1,1)
        imagesc(stackedDataSub(:,:,ind));
        colorbar()
        drawnow
        
        subplot(3,1,2)
        contour(stackedDataSub(:,:,ind));
        set(gca,'YDir','reverse')
        colorbar()
        drawnow()
        
        subplot(3,1,3)
        histogram(reshape(stackedDataSub(:,:,ind),[],1),100,'Normalization','probability')
        drawnow
        
    end
end
%% 3D plot to look at everything
if plotIt
    
    fig3 = figure;
    for ind = 1:size(stackedDataSub,3)
        surf(stackedDataSub(:,:,ind));
        switch file
            case 1
                zlim([23.5 25.5])
            case 4
                zlim([22.4 23.8])
            case 5
                zlim([21.4 22.8])
            case 6
                zlim([22.1 23])
        end
        drawnow
    end
    
end
%% gradient

[g_x,g_y] = gradient(stackedDataSub);

if plotIt
    
    fig4 = figure;
    for ind = 1:size(stackedDataSub,3)
        
        quiver(g_x(:,:,ind),g_y(:,:,ind))
        set(gca,'YDir','reverse')
        xlim([1 121])
        ylim([1 171])
        
        drawnow
        
    end
    
end

%% find max

[val,ind] = max(stackedDataSub(:));
[x_max,y_max,t_max] = ind2sub(size(stackedDataSub),ind);
%% visualize
trialInt = t_max;
dataInt = stackedDataSub(:,:,trialInt);
figure
contour(dataInt)
set(gca,'YDir','reverse')

figure
surf(dataInt);
set(gca,'YDir','reverse')


%% experimental analysis
% SVD
unwrapped = reshape(stackedDataSub,[],size(stackedDataSub,3));

% get index of interest for to check reshaping
% switch file
%     case 3
%         pix_int = sub2ind(size(stackedDataSub(:,:,1)),71,48); % case 3
%     case 4
%         pix_int = sub2ind(size(stackedDataSub(:,:,1)),94,92); % case 4 % peak trial 1172
% end

pix_int = sub2ind(size(stackedDataSub(:,:,1)),x_max,y_max);
figure
plot(unwrapped(pix_int,:))

%%
[U,S,V] = svd(unwrapped','econ');
%%
figure
plot((diag(S)/sum(diag(S))),'o')

figure
for i = 1:3
    subplot(3,1,i)
    plot(U(:,i))
end

figure
for i = 1:3
    subplot(3,1,i)
    plot(V(:,i))
end

figure
for i = 1:3
    subplot(3,1,i)
    imagesc(reshape(V(:,i),size(stackedDataSub,1),size(stackedDataSub,2)))
    
end
%%
% 3 mode reconstruction
begin_mode = 2;
end_mode = 4;

recon = U(:,begin_mode:end_mode)*S(begin_mode:end_mode,begin_mode:end_mode)*V(:,begin_mode:end_mode)';
recon = recon';
rewrapped = reshape(recon,size(stackedDataSub,1),size(stackedDataSub,2),size(stackedDataSub,3));

if plotIt
    
    figure;
    
    for ind = 1:size(rewrapped,3)
        
        subplot(3,1,1)
        imagesc(rewrapped(:,:,ind));
        colorbar()
        drawnow
        
        subplot(3,1,2)
        contour(rewrapped(:,:,ind));
        set(gca,'YDir','reverse')
        colorbar()
        drawnow()
        
        subplot(3,1,3)
        hist(reshape(rewrapped(:,:,ind),[],1),100)
        drawnow
        
    end
end

%% 3D plot to look at everything
if plotIt
    
    figure;
    for ind = 1:size(rewrapped,3)
        surf(rewrapped(:,:,ind));
        drawnow
    end
    
end

%% - DMD
[Phi, mu, lambda, diagS, x0] = DMD(unwrapped);
