% 7-14-2017 - script to analyze heating profile in the thermal camera
% imaging data

%% load file
close all;clear all;clc

file = 1;
switch file
    case 1
        xVals = [240:380];
        yVals = [180:300];
        load('C:\david\david\11-8-2017\2nd_gel_flipped_betweenGrids\2ndGel_betweenGrids_29_32_run1_converted\2ndGel_betweenGrids_29_32_run1.mat')
        stackedData = stackedData - 273.15;
        
        stackedDataSub1 = stackedData(xVals,yVals,:);
        stackedDataSub1 = baselineSubtract(stackedDataSub1);
        
        clear stackedData;
        load('C:\david\david\11-8-2017\2nd_gel_flipped_betweenGrids\2ndGel_betweenGrids_29_32_run2_converted\2ndGel_betweenGrids_29_32_run2.mat')
        stackedData = stackedData - 273.15;
        
        stackedDataSub2 = stackedData(xVals,yVals,:);
        stackedDataSub2 = baselineSubtract(stackedDataSub2);
        
        clear stackedData;
        load('C:\david\david\11-8-2017\2nd_gel_flipped_betweenGrids\2ndGel_betweenGrids_29_32_run3_converted\2ndGel_betweenGrids_29_32_run3.mat')
        stackedData = stackedData - 273.15;
        
        stackedDataSub3 = stackedData(xVals,yVals,:);
        stackedDataSub3 = baselineSubtract(stackedDataSub3);
        
        clear stackedData;
        
        stackedDataSub = cat(4,stackedDataSub1(:,:,1:end-1), stackedDataSub2(:,:,1:end-1), stackedDataSub3);
        %stackedDataSub = cat(4,stackedDataSub1(:,:,1:end-1), stackedDataSub2(:,:,1:end-1));

    case 2
        xVals = [200:400];
        yVals = [200:400];
        load('C:\david\david\10-8-2017\optistim4_dual_1.75_3.05_11_19_27_trial1_converted\optistim4_dual_1.75_3.05_11_19_27_trial1.mat')
        stackedDataSub1 = stackedData(xVals,yVals,:);
        stackedDataSub1 = baselineSubtract(stackedDataSub1);
        
        clear stackedData;
        load('C:\david\david\10-8-2017\optistim4_dual_1.75_3.05_11_19_27_trial2_converted\optistim4_dual_1.75_3.05_11_19_27_trial2.mat')
        
        stackedDataSub2 = stackedData(xVals,yVals,:);
        % mean subtract first frame
        stackedDataSub2 = baselineSubtract(stackedDataSub2);
        clear stackedData;
        load('C:\david\david\10-8-2017\optistim4_dual_1.75_3.05_11_19_27_trial3_converted\optistim4_dual_1.75_3.05_11_19_27_trial3.mat')
        stackedDataSub3 = stackedData(xVals,yVals,:);
        stackedDataSub3 = baselineSubtract(stackedDataSub3);
        clear stackedData;
        
        
        stackedDataSub = cat(4,stackedDataSub1, stackedDataSub2, stackedDataSub3);
        
    case 3
        xVals = [200:400];
        yVals = [200:400];
        load('C:\david\david\10-8-2017\optistim8_long_1.75_3.05_4_29_1_converted\optstim8_long_1.75_3.05_4_29_1.mat')
        stackedDataSub1 = stackedData(xVals,yVals,:);
        stackedDataSub1 = baselineSubtract(stackedDataSub1);
        
        clear stackedData;
        load('C:\david\david\10-8-2017\optistim8_long_1.75_3.05_4_29_2_converted\optstim8_long_1.75_3.05_4_29_2.mat')
        
        stackedDataSub2 = stackedData(xVals,yVals,:);
        % mean subtract first frame
        stackedDataSub2 = baselineSubtract(stackedDataSub2);
        clear stackedData;
        
        stackedDataSub = cat(4,stackedDataSub1, stackedDataSub2);
    case 4
        xVals = [200:400];
        yVals = [1:150];
        load('C:\david\david\10-8-2017\optstim8_1.75_3.05_4_29_1_converted\optstim8_1.75_3.05_4_29_1.mat')
        stackedDataSub1 = stackedData(xVals,yVals,:);
        stackedDataSub1 = baselineSubtract(stackedDataSub1);
        
        clear stackedData;
        load('C:\david\david\10-8-2017\optstim8_1.75_3.05_4_29_2_converted\optstim8_1.75_3.05_4_29_2.mat')
        
        stackedDataSub2 = stackedData(xVals,yVals,:);
        % mean subtract first frame
        stackedDataSub2 = baselineSubtract(stackedDataSub2);
        clear stackedData;
        load('C:\david\david\10-8-2017\optstim8_1.75_3.05_4_29_3_converted\optstim8_1.75_3.05_4_29_3.mat')
        stackedDataSub3 = stackedData(xVals,yVals,:);
        stackedDataSub3 = baselineSubtract(stackedDataSub3);
        clear stackedData;
        
        
        stackedDataSub = cat(4,stackedDataSub1, stackedDataSub2, stackedDataSub3(:,:,1:end-1));
end
%%
stackedDataSubAverage = mean(stackedDataSub,4);

smoothIt = 1;

% 2d smoothing
if smoothIt
    %h = 1/9*[0 1 0; 1 5 1; 0 1 0];
    h = 0.125*ones(3);
    for i = 1:size(stackedDataSubAverage,3)
        stackedDataSubAverage(:,:,i) = conv2(stackedDataSubAverage(:,:,i),h,'same');
    end
    
    stackedDataSubAverage = stackedDataSubAverage(2:end-1,2:end-1,:);
end

plotIt = 1;

%% colorbrewer palatte

CT=cbrewer('seq', 'Blues', 8);
colormap(CT)
%% begin the analysis

% visualize
if plotIt
    
    fig1 = figure;
    for ind = 1:size(stackedData,3)
        
        subplot(3,1,1)
        imagesc(stackedDataAverage(:,:,ind));
        colorbar()
        drawnow()
        
        subplot(3,1,2)
        contour(stackedDataAverage(:,:,ind));
        set(gca,'YDir','reverse')
        colorbar()
        drawnow()
        
        subplot(3,1,3)
        histogram(reshape(stackedDataAverage(:,:,ind),[],1),100,'Normalization','probability')
        drawnow()
        
    end
end
%% look at subportion of thermal heating



if plotIt
    
    fig2 = figure;
    
    for ind = 1:size(stackedDataSubAverage,3)
        
        subplot(3,1,1)
        imagesc(stackedDataSubAverage(:,:,ind));
        colorbar()
        drawnow
        
        subplot(3,1,2)
        contour(stackedDataSubAverage(:,:,ind));
        set(gca,'YDir','reverse')
        colorbar()
        drawnow()
        
        subplot(3,1,3)
        histogram(reshape(stackedDataSubAverage(:,:,ind),[],1),100,'Normalization','probability')
        drawnow
        
    end
end
%% 3D plot to look at everything
if plotIt
    
    fig3 = figure;
    CT=cbrewer('seq', 'Reds', 20);
    colormap(CT)
    for ind = 1:10:size(stackedDataSubAverage,3)
        surf(stackedDataSubAverage(:,:,ind));
        
        switch file
            case 1
                %                 zlim([22.1 23])
                zlim([-0.25 0.75])
                view([-77.10, 57.2])
                
            case 2
                zlim([-1 2])
            case 3
                zlim([-0.25 0.9])
                view([31.7,62.8])
            case 4
                zlim([-0.25 0.9])
                
        end
        drawnow
    end
    
end

%%
if plotIt
    fig3 = figure;
    CT=cbrewer('seq', 'Reds', 20);
    colormap(CT)
    for ind = 1:10:size(stackedDataSubAverage,3)
        imagesc(stackedDataSubAverage(:,:,ind));
        
        colorbar()
        
        switch file
            case 1
                caxis([-0.25 0.75])
                caxis([0 0.75])
                
            case 2
                zlim([-1 2])
            case 3
                caxis([-0.25 0.9])
            case 4
                caxis([-0.25 0.9])
                
        end
        drawnow
    end
    
end

%% example case of
ind = 2161; % works well for elliptical, 
figure
zlim([22.1 23])

subplot(2,1,1)
surf(stackedDataSub1(:,:,ind))
subplot(2,1,2)
surf(stackedDataSubAverage(:,:,ind))
subplot(2,1,1)
title('Individual Recording')
zlim([-0.2 0.6])
subplot(2,1,2)
title('Average Recording')
zlabel('Temperature (Celsius)')
zlim([-0.2 0.6])


%% compare noise levels
figure
histogram(stackedDataSub1(:,:,ind),'normalization','probability')
hold on
histogram(stackedDataSubAverage(:,:,ind),'normalization','probability')
legend({'Individual Trial','Average'})
title('Distribution of temperature values in frame')
%% gradient

[g_x,g_y] = gradient(stackedDataSubAverage);

if plotIt
    
    fig4 = figure;
    for ind = 1:size(stackedDataSubAverage,3)
        
        quiver(g_x(:,:,ind),g_y(:,:,ind))
        set(gca,'YDir','reverse')
        xlim([1 121])
        ylim([1 171])
        
        drawnow
        
    end
    
end

%% find max

[val,ind] = max(stackedDataSubAverage(:));
[x_max,y_max,t_max] = ind2sub(size(stackedDataSubAverage),ind);


%% visualize
trialInt = t_max;
dataInt = stackedDataSubAverage(:,:,trialInt);
figure
imagesc(dataInt)
colormap('hot')
set(gca,'YDir','reverse')
h = colorbar;
ylabel(h, ['Temperature Change in ' char(176) 'K'])
ylabel(['Pixel'])
xlabel(['Pixel'])
title(['Temperature difference from baseline at maximum heating point ' char(176) 'K'])

figure
surf(dataInt);
colormap('hot')

set(gca,'YDir','reverse')
zlabel(['Temperature Change in ' char(176) 'K']) 
ylabel(['Pixel'])
xlabel(['Pixel'])
title(['Temperature difference from baseline at maximum heating point ' char(176) 'K'])

%% fit gaussian

% separate peaks
switch file
    case 1
        center1 = [68 60];
        center2 = [62 95];
        pixels = 20;
        area1 = dataInt(center1(2)-2*pixels:center1(2)+pixels/2,center1(1)-1.5*pixels:center1(1)+1.5*pixels);
area2 = dataInt(center2(2)-pixels/2:center2(2)+2*pixels,center2(1)-1.5*pixels:center2(1)+1.5*pixels);
        % calculate distance between centers
        distance = norm(center2 - center1);
        pixel_dist_conversion = 0.2640;
        
    case 2
        center1 = [86 119];
        pixels = 40;
                area1 = dataInt(center1(2)-pixels:center1(2)+pixels,center1(1)-pixels:center1(1)+pixels);
%area2 = dataInt(center2(2)-pixels:center2(2)+pixels,center2(1)-pixels:center2(1)+pixels);

end
figure
imagesc(area1)
fit1 = D2gaussianFit_function(area1);

%%
figure
imagesc(area2)
fit2 = D2gaussianFit_function(area2);

% time to compare the fits of these

x0 = 1;
y0 = 1;
y1 = size(dataInt,1);
x1 = size(dataInt,2);
grid_data = makeGrid(x0,x1,y0,y1);
% recon 1
f_out1 = D2GaussFunctionRot_arbCenter(fit1,grid_data,center1(1),center1(2));

% recon 2
f_out2 = D2GaussFunctionRot_arbCenter(fit2,grid_data,center2(1),center2(2));

figure
C = del2(dataInt);
surf(dataInt,C)
hold on
surf(f_out1);
alpha(0.4)
colormap('hot')
hold on
surf(f_out2);
alpha(0.4)
colormap('hot')
zlabel('Temperature difference from baseline (Kelvin)')


%% figure out time lag

for i = 1:size(stackedDataSub,4)
    unwrapped(:,:,i) = reshape(stackedDataSub(:,:,:,i),[],size(stackedDataSub,3));
end

% get index of interest for to check reshaping
% switch file
%     case 3
%         pix_int = sub2ind(size(stackedDataSub(:,:,1)),71,48); % case 3
%     case 4
%         pix_int = sub2ind(size(stackedDataSub(:,:,1)),94,92); % case 4 % peak trial 1172
% end

for i = 1:size(unwrapped,3)
    [val,ind] = max(reshape(stackedDataSub(:,:,:,i),1,[]));
    [x_max(i),y_max(i),t_max(i)] = ind2sub(size(stackedDataSub(:,:,:,i)),ind);
    pix_int(:,i) = sub2ind(size(stackedDataSub(:,:,1,i)),x_max(i),y_max(i));
    pixel(:,i) = unwrapped(pix_int(:,i),:,i);
end


figure
for i = 1:size(pixel,2)
    plot(pixel(:,i))
    hold on
end

%% try smoothing it
%
order = 3;
framelen = 321;
framelen = 721;

sgf = sgolayfilt(pixel,order,framelen);

close all;
figure
for i = 1:size(pixel,2)
    plot(pixel(:,i))
    hold on
    plot(sgf(:,i),'linewidth',2,'color','k')
    legend('signal','sgolay')
    ylim([-0.2 1])
    xlabel('sample number')
    ylabel(['temperature change from baseline ',char(176),'K'])
end
title('Savitsky-Golay filtered time courses of temperature profiles')
%% xcorr

%sgf
Fs = 30;
[cr,lgs] = xcorr(sgf,'coeff');
numSigs = size(sgf,2);

figure
for row = 1:numSigs;
    for col = 1:numSigs
        nm = numSigs*(row-1)+col;
        subplot(numSigs,numSigs,nm)
        plot(lgs,cr(:,nm),'.')
        title(sprintf('c_{%d%d}',row,col))
        xlim([-500 500])
    end
    
end
    xlabel('sample number')
    ylabel('normalized cross-correlation')

% align relative to the first one 

[~,I] = max(abs(cr));
lagDiff = lgs(I);
fprintf('the lags are %d samples \n ',lagDiff(1:numSigs))
timeDiff = lagDiff/Fs;
fprintf('the lags are %f seconds \n ',timeDiff(1:numSigs))

%% more experimental analysis. compare FEM results to the data, see the fit

[x,y,T] = importfile_quickfield('simpleHeating_rightSize_smallerCurrent.txt');

%%
% trim down simulation 
T = T(:,1:105);
x = x(:,1:105);
y = y(:,1:105);

sizeM_x = size(x,1);
sizeM_y = size(x,2);

dataInt_s = dataInt'+20;
dataInt_s = dataInt_s(35:35+sizeM_x-1,15:15+sizeM_y-1);
[d,Z,transform] = procrustes(T,dataInt_s);
%%
figure
subplot(3,1,1)
imagesc(T)
colorbar()
title('FEM model')

subplot(3,1,2)
imagesc(dataInt_s)
colorbar()
title('Original Data')

subplot(3,1,3)
imagesc(Z)
title('Procrustes transformed data')
colorbar()
%%
subplot(2,1,1);
imagesc(dataInt_s)
subplot(2,1,2)
imagesc(T)
%%
Fs = 30;
[acor,lag] = xcorr(pix_int1,pix_int2);

[~,I] = max(abs(acor));
lagDiff = lag(I);
timeDiff = lagDiff/Fs;

figure
plot(lag,acor)

Fs = 30;
[acor,lag] = xcorr(pix_int1,pix_int2);

[~,I] = max(abs(acor));
lagDiff = lag(I);
timeDiff = lagDiff/Fs;

figure
plot(lag,acor)


%%
s1al = s1(-lagDiff:end);
t1al = (0:length(s1al)-1)/Fs;

subplot(2,1,1)
plot(t1al,s1al)
title('s_1, aligned')

subplot(2,1,2)
plot(t2,s2)
title('s_2')
xlabel('Time (s)')


%% experimental analysis
% SVD
unwrapped = reshape(stackedDataSubAverage,[],size(stackedDataSubAverage,3));

% get index of interest for to check reshaping
% switch file
%     case 3
%         pix_int = sub2ind(size(stackedDataSub(:,:,1)),71,48); % case 3
%     case 4
%         pix_int = sub2ind(size(stackedDataSub(:,:,1)),94,92); % case 4 % peak trial 1172
% end

pix_int = sub2ind(size(stackedDataSubAverage(:,:,1)),x_max,y_max);
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
    imagesc(reshape(V(:,i),size(stackedDataSubAverage,1),size(stackedDataSubAverage,2)))
    
end
%%
% 3 mode reconstruction
begin_mode = 2;
end_mode = 4;

recon = U(:,begin_mode:end_mode)*S(begin_mode:end_mode,begin_mode:end_mode)*V(:,begin_mode:end_mode)';
recon = recon';
rewrapped = reshape(recon,size(stackedDataSubAverage,1),size(stackedDataSubAverage,2),size(stackedDataSubAverage,3));

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
