% The prototype for interactive matrix profile calculation
% Chin-Chia Michael Yeh 01/26/2016
%
% [matrixProfile, profileIndex, motifIndex] = interactiveMatrixProfile(data, subsequenceLength);
%

function [matrixProfile, profileIndex, motifIdx] = interactiveMatrixProfile(data, subLen)
%% set trivial match exclusion zone
exclusionZone = round(subLen/2);
% exclusionZone = round(subLen/4);

%% check input
dataLen = length(data);
if subLen > dataLen/2
    error('Error: Time series is too short relative to desired subsequence length');
end
if subLen < 4
    error('Error: Subsequence length must be at least 4');
end
if dataLen == size(data, 2)
    data = data';
end

%% spawn main window
mainWindow.fig = figure('name', 'Interactive Matrix Profile Calculation', ...
    'visible', 'off', 'toolbar', 'none', 'ResizeFcn', @mainResize);

%% add UI element into the window
figPosition = get(mainWindow.fig, 'position');
backColor = get(mainWindow.fig, 'color');
axesHeight = round((figPosition(4)-140)/3);
mainWindow.dataAx = axes('parent',mainWindow.fig, 'units', 'pixels', ...
    'position', [30, 2*axesHeight+110, figPosition(3)-60, axesHeight]);
mainWindow.profileAx = axes('parent',mainWindow.fig, 'units', 'pixels', ...
    'position', [30, axesHeight+70, figPosition(3)-60, axesHeight]);
mainWindow.motifAx = axes('parent',mainWindow.fig, 'units', 'pixels', ...
    'position', [30, 30, figPosition(3)-160, axesHeight]);
mainWindow.discardBtn = uicontrol('parent',mainWindow.fig, 'style', 'pushbutton',...
    'string', 'Discard', 'fontsize', 10, 'position', [figPosition(3)-120, 65, 90, 25], ...
    'callback', @pushDiscardBtn);
mainWindow.stopBtn = uicontrol('parent',mainWindow.fig, 'style', 'pushbutton',...
    'string', 'Stop', 'fontsize', 10, 'position', [figPosition(3)-120, 30, 90, 25], ...
    'callback', @pushStopBtn);
mainWindow.dataText = uicontrol('parent',mainWindow.fig, 'style', 'text',...
    'string', '', 'fontsize', 10, 'backgroundcolor', backColor, ...
    'horizontalalignment', 'left', 'position', [30, 3*axesHeight+110, figPosition(3)-60, 18]);
mainWindow.motifText = uicontrol('parent',mainWindow.fig, 'style', 'text',...
    'string', '', 'fontsize', 10, 'backgroundcolor', backColor, ...
    'horizontalalignment', 'left', 'position', [30, axesHeight+30, figPosition(3)-160, 18]);
mainWindow.profileText = uicontrol('parent',mainWindow.fig, 'style', 'text',...
    'string', 'The best-so-far matrix profile', 'fontsize', 10, ...
    'backgroundcolor', backColor,'horizontalalignment', 'left', ...
    'position', [30, 2*axesHeight+70, figPosition(3)-60, 18]);

%% modify the properties of the axis
set(mainWindow.dataAx,'xlim',[1, dataLen]);
set(mainWindow.dataAx,'ylim',[-0.05, 1.05]);
set(mainWindow.dataAx,'ytick',[]);
set(mainWindow.dataAx,'ycolor',[1 1 1]);
set(mainWindow.profileAx,'xlim',[1, dataLen]);
set(mainWindow.profileAx,'ylim',[0, 2*sqrt(subLen)]);
set(mainWindow.motifAx,'xlim',[1, subLen]);
set(mainWindow.motifAx,'ylim',[-0.05, 1.05]);
set(mainWindow.motifAx,'ytick',[]);
set(mainWindow.motifAx,'ycolor',[1 1 1]);

%% plot data
dataPlot = zeroOneNorm(data);
hold(mainWindow.dataAx, 'on');
plot(1:dataLen, dataPlot, 'r', 'parent', mainWindow.dataAx);
hold(mainWindow.dataAx, 'off');

%% preprocess for matrix profile
[dataFreq, data2Sum, dataSum, dataMean, data2Sig, dataSig] = ...
    fastfindNNPre(data, dataLen, subLen);
profileLen = dataLen - subLen + 1;
idxOrder = randperm(profileLen);
matrixProfile = inf(profileLen, 1);
profileIndex = zeros(profileLen, 1);

%% iteratively plot
mainWindow.stopping = false;
mainWindow.discardIdx = [];
set(mainWindow.fig, 'userdata', mainWindow);
firstUpdate = true;
timer = tic();
for i = 1:profileLen
    % compute the distance profile
    idx = idxOrder(i);
    query = data(idx:idx+subLen-1);
    distanceProfile = fastfindNN(dataFreq, query, dataLen, subLen, ...
        data2Sum, dataSum, dataMean, data2Sig, dataSig);
    distanceProfile = abs(distanceProfile);
    
    % apply exclusion zone
    exclusionZoneStart = max(1, idx-exclusionZone);
    exclusionZoneEnd = min(profileLen, idx+exclusionZone);
    distanceProfile(exclusionZoneStart:exclusionZoneEnd) = inf;
    
    % figure out and store the neareest neighbor
    if i == 1
        matrixProfile = distanceProfile;
        profileIndex(:) = idx;
    else
        updatePos = distanceProfile < matrixProfile;
        profileIndex(updatePos) = idx;
        matrixProfile(updatePos) = distanceProfile(updatePos);
    end
    [matrixProfile(idx), profileIndex(idx)] = min(distanceProfile);
    
    % plotting
    if toc(timer) > 1 || i == profileLen
        % plot matrix profile
        if exist('prefilePlot', 'var')
            delete(prefilePlot);
        end
        hold(mainWindow.profileAx, 'on');
        prefilePlot = plot(1:profileLen, matrixProfile, 'b', 'parent', mainWindow.profileAx);
        hold(mainWindow.profileAx, 'off');
        
        % plot motif
        if exist('motifDataPlot', 'var')
            for j = 1:2
                delete(motifDataPlot(j));
            end
        end
        if exist('motifMotifPlot', 'var')
            for j = 1:2
                delete(motifMotifPlot(j));
            end
        end
        mainWindow = get(mainWindow.fig, 'userdata');
        discardIdx = mainWindow.discardIdx;
        matrixProfileTemp = matrixProfile;
        for j = 1:length(discardIdx)
            discardZoneStart = max(1, discardIdx(j)-exclusionZone);
            discardZoneEnd = min(profileLen, discardIdx(j)+exclusionZone);
            matrixProfileTemp(discardZoneStart:discardZoneEnd) = inf;
            matrixProfileTemp(abs(profileIndex - discardIdx(j)) < exclusionZone) = inf;
        end
        [~, minIdx] = min(matrixProfileTemp);
        motifIdx = sort([minIdx, profileIndex(minIdx)]);
        motifDataPlot = zeros(2, 1);
        motifMotifPlot = zeros(2, 1);
        motifColor = {'g', 'c'};
        for j = 1:2
            motifPos = motifIdx(j):motifIdx(j)+subLen-1;
            hold(mainWindow.dataAx, 'on');
            motifDataPlot(j) = plot(motifPos, dataPlot(motifPos), ...
                motifColor{j}, 'parent', mainWindow.dataAx);
            hold(mainWindow.dataAx, 'off');
            hold(mainWindow.motifAx, 'on');
            motifMotifPlot(j) = plot(1:subLen, zeroOneNorm(dataPlot(motifPos)),...
                motifColor{j}, 'parent', mainWindow.motifAx);
            hold(mainWindow.motifAx, 'off');
        end
        set(mainWindow.dataText, 'string', ...
            sprintf('We are %.1f%% done: The input time series: The best-so-far motifs are color coded (see bottom panel)', i*100/profileLen));
        set(mainWindow.motifText, 'string', ...
            sprintf('The best-so-far motifs are located at %d (green) and %d (cyan)', motifIdx(1), motifIdx(2)));
        
        % show the figure
        if firstUpdate
            set(mainWindow.fig, 'userdata', mainWindow);
            set(mainWindow.fig, 'visible', 'on');
            firstUpdate = false;
        end
        
        % check for stop
        mainWindow = get(mainWindow.fig, 'userdata');
        mainWindow.motifIdx = motifIdx;
        set(mainWindow.fig, 'userdata', mainWindow);
        if mainWindow.stopping
            set(mainWindow.fig, 'name', 'Interactive Matrix Profile Calculation (Stopped)');
            return;
        end
        if i == profileLen
            set(mainWindow.fig, 'name', 'Interactive Matrix Profile Calculation (Completed)');
            set(mainWindow.discardBtn, 'enable', 'off');
            set(mainWindow.stopBtn, 'enable', 'off');
            return;
        end
        
        % pause for plot and restart timer
        pause(eps);
        timer = tic();
    end
end

function pushDiscardBtn(src, ~)
mainWindowFig = get(src, 'parent');
mainWindow = get(mainWindowFig, 'userdata');
mainWindow.discardIdx = [mainWindow.discardIdx, mainWindow.motifIdx];
set(mainWindow.fig, 'userdata', mainWindow);

function pushStopBtn(src, ~)
mainWindowFig = get(src, 'parent');
mainWindow = get(mainWindowFig, 'userdata');
mainWindow.stopping = true;
set(mainWindow.discardBtn, 'enable', 'off');
set(src, 'enable', 'off');
set(mainWindow.fig, 'userdata', mainWindow);

function [dataFreq, data2Sum, dataSum, dataMean, data2Sig, dataSig] = ...
    fastfindNNPre(data, dataLen, subLen)
data(dataLen+1:2*dataLen) = 0;
dataFreq = fft(data);
cum_sumx = cumsum(data);
cum_sumx2 =  cumsum(data.^2);
data2Sum = cum_sumx2(subLen:dataLen)-[0;cum_sumx2(1:dataLen-subLen)];
dataSum = cum_sumx(subLen:dataLen)-[0;cum_sumx(1:dataLen-subLen)];
dataMean = dataSum./subLen;
data2Sig = (data2Sum./subLen)-(dataMean.^2);
dataSig = sqrt(data2Sig);

function distanceProfile = fastfindNN(dataFreq, query, dataLen, subLen, ...
    data2Sum, dataSum, dataMean, data2Sig, dataSig)
query = (query-mean(query))./std(query,1);
query = query(end:-1:1);
query(subLen+1:2*dataLen) = 0;
queryFreq = fft(query);
dataQueryProdFreq = dataFreq.*queryFreq;
dataQueryProd = ifft(dataQueryProdFreq);
querySum = sum(query);
query2Sum = sum(query.^2);
distanceProfile = (data2Sum - 2*dataSum.*dataMean + subLen*(dataMean.^2))./data2Sig ...
    - 2*(dataQueryProd(subLen:dataLen) - querySum.*dataMean)./dataSig + query2Sum;
distanceProfile = sqrt(distanceProfile);

function x = zeroOneNorm(x)
x = x-min(x);
x = x/max(x);

function mainResize(src, ~)
mainWindow = get(src, 'userdata');
figPosition = get(mainWindow.fig, 'position');
axesHeight = round((figPosition(4)-140)/3);
set(mainWindow.dataAx, 'position', [30, 2*axesHeight+110, figPosition(3)-60, axesHeight]);
set(mainWindow.profileAx, 'position', [30, axesHeight+70, figPosition(3)-60, axesHeight]);
set(mainWindow.motifAx, 'position', [30, 30, figPosition(3)-160, axesHeight]);
set(mainWindow.discardBtn, 'position', [figPosition(3)-120, 65, 90, 25]);
set(mainWindow.stopBtn, 'position', [figPosition(3)-120, 30, 90, 25]);
set(mainWindow.dataText, 'position', [30, 3*axesHeight+110, figPosition(3)-60, 18]);
set(mainWindow.profileText, 'position', [30, 2*axesHeight+70, figPosition(3)-60, 18]);
set(mainWindow.motifText, 'position', [30, axesHeight+30, figPosition(3)-160, 18]);