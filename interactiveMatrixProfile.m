% The prototype for interactive matrix profile calculation
% Chin-Chia Michael Yeh 01/26/2016
%
% [matrixProfile, profileIndex, motifIndex, discordIndex] = ...
%     interactiveMatrixProfile(data, subLen);
% Output:
%     matrixProfile: matrix porfile of the self-join (vector)
%     profileIndex: matrix porfile index of the self-join (vector)
%     motifIndex: index of the first, second, and third motifs and their associated nearest neighbors when stopped (3x2 cell)
%                +--------------------------------+-------------------------------------------+
%                | pair of index for first motif  | nearest neighbor of the first motif pair  |
%                +--------------------------------+-------------------------------------------+
%                | pair of index for second motif | nearest neighbor of the second motif pair |
%                +--------------------------------+-------------------------------------------+
%                | pair of index for third motif  | nearest neighbor of the third motif pair  |
%                +--------------------------------+-------------------------------------------+
%     discordIndex: index of discords when stopped (vector)
% Input:
%     data: input time series (vector)
%     subLen: interested subsequence length (scalar)
%
% Chin-Chia Michael Yeh, Yan Zhu, Liudmila Ulanova, Nurjahan Begum, Yifei Ding, Hoang Anh Dau, 
% Diego Furtado Silva, Abdullah Mueen, and Eamonn Keogh, "Matrix Profile I: All Pairs Similarity 
% Joins for Time Series," ICDM 2016, http://www.cs.ucr.edu/~eamonn/MatrixProfile.html
%

function [matrixProfile, profileIndex, motifIdxs, discordIdx] = ...
    interactiveMatrixProfile(data, subLen)
%% set trivial match exclusion zone and motif radius
exclusionZone = round(subLen * 0.5);
radius = 2;

%% check input
dataLen = length(data);
if subLen > dataLen / 2
    error(['Error: Time series is too short ', ...
        'relative to desired subsequence length']);
end
if subLen < 4
    error('Error: Subsequence length must be at least 4');
end
if subLen > dataLen / 20
    error('Error: subsequenceLength > dataLength / 20')
end
if dataLen == size(data, 2)
    data = data';
end

%% spawn main window
mainWindow.fig = figure('name', ...
    'UCR Interactive Matrix Profile Calculation', ...
    'visible', 'off', 'toolbar', 'none', 'ResizeFcn', @mainResize);

%% add UI element into the window
backColor = get(mainWindow.fig, 'color');
mainWindow.dataAx = axes('parent', mainWindow.fig, 'units', 'pixels');
mainWindow.profileAx = axes('parent', mainWindow.fig, 'units', 'pixels');
mainWindow.motif1Ax = axes('parent', mainWindow.fig, 'units', 'pixels');
mainWindow.motif2Ax = axes('parent', mainWindow.fig, 'units', 'pixels');
mainWindow.motif3Ax = axes('parent', mainWindow.fig, 'units', 'pixels');
mainWindow.discordAx = axes('parent', mainWindow.fig, 'units', 'pixels');
mainWindow.discardBtn = zeros(3, 1);
mainWindow.discard1Btn = uicontrol('parent',mainWindow.fig, ...
    'style', 'pushbutton', 'string', 'Discard', 'fontsize', 10, ...
    'callback', @(src, cbdata) pushDiscardBtn(src, cbdata, 1));
mainWindow.discard2Btn = uicontrol('parent', mainWindow.fig, ...
    'style', 'pushbutton', 'string', 'Discard', 'fontsize', 10, ...
    'callback', @(src, cbdata) pushDiscardBtn(src, cbdata, 2));
mainWindow.discard3Btn = uicontrol('parent', mainWindow.fig, ...
    'style', 'pushbutton', 'string', 'Discard', 'fontsize', 10, ...
    'callback', @(src, cbdata) pushDiscardBtn(src, cbdata, 3));
mainWindow.stopBtn = uicontrol('parent', mainWindow.fig, ...
    'style', 'pushbutton', 'string', 'Stop', 'fontsize', 10, ...
    'callback', @pushStopBtn);
mainWindow.dataText = uicontrol('parent', mainWindow.fig, ...
    'style', 'text', 'string', '', 'fontsize', 10, ...
    'backgroundcolor', backColor, 'horizontalalignment', 'left');
mainWindow.profileText = uicontrol('parent', mainWindow.fig, ...
    'style', 'text', 'string', 'The best-so-far matrix profile', ...
    'fontsize', 10, 'backgroundcolor', backColor, ...
    'horizontalalignment', 'left');
mainWindow.motif1Text = uicontrol('parent', mainWindow.fig, ...
    'style', 'text', 'string', '', 'fontsize', 10, ...
    'backgroundcolor', backColor, 'horizontalalignment', 'left');
mainWindow.motif2Text = uicontrol('parent', mainWindow.fig, ...
    'style', 'text', 'string', '', 'fontsize', 10, ...
    'backgroundcolor', backColor, 'horizontalalignment', 'left');
mainWindow.motif3Text = uicontrol('parent', mainWindow.fig, ...
    'style', 'text', 'string', '', 'fontsize', 10, ...
    'backgroundcolor', backColor, 'horizontalalignment', 'left');
mainWindow.discordText = uicontrol('parent', mainWindow.fig, ...
    'style', 'text', 'string', '', 'fontsize', 10, ...
    'backgroundcolor', backColor, 'horizontalalignment', 'left');

%% modify the properties of the axis
set(mainWindow.dataAx,'xlim', [1, dataLen]);
set(mainWindow.dataAx,'xtick', [1, dataLen]);
set(mainWindow.dataAx,'ylim', [-0.05, 1.05]);
set(mainWindow.dataAx,'ytick', []);
set(mainWindow.dataAx,'ycolor', [1 1 1]);
set(mainWindow.profileAx,'xlim', [1, dataLen]);
set(mainWindow.profileAx,'xtick', [1, dataLen]);
set(mainWindow.profileAx,'ylim', [0, 2*sqrt(subLen)]);
set(mainWindow.motif1Ax,'xlim', [1, subLen]);
set(mainWindow.motif1Ax,'xtick', [1, subLen]);
set(mainWindow.motif1Ax,'ylim', [-0.05, 1.05]);
set(mainWindow.motif1Ax,'ytick', []);
set(mainWindow.motif1Ax,'ycolor', [1 1 1]);
set(mainWindow.motif2Ax,'xlim', [1, subLen]);
set(mainWindow.motif2Ax,'xtick', [1, subLen]);
set(mainWindow.motif2Ax,'ylim', [-0.05, 1.05]);
set(mainWindow.motif2Ax,'ytick', []);
set(mainWindow.motif2Ax,'ycolor', [1 1 1]);
set(mainWindow.motif3Ax,'xlim', [1, subLen]);
set(mainWindow.motif3Ax,'xtick', [1, subLen]);
set(mainWindow.motif3Ax,'ylim', [-0.05, 1.05]);
set(mainWindow.motif3Ax,'ytick', []);
set(mainWindow.motif3Ax,'ycolor', [1 1 1]);
set(mainWindow.discordAx,'xlim', [1, subLen]);
set(mainWindow.discordAx,'xtick', [1, subLen]);
set(mainWindow.discordAx,'ylim', [-0.05, 1.05]);
set(mainWindow.discordAx,'ytick', []);
set(mainWindow.discordAx,'ycolor', [1 1 1]);

%% plot data
dataPlot = zeroOneNorm(data);
hold(mainWindow.dataAx, 'on');
plot(1:dataLen, dataPlot, 'r', 'parent', mainWindow.dataAx);
hold(mainWindow.dataAx, 'off');

%% locate nan and inf
profileLen = dataLen - subLen + 1;
isSkip = false(profileLen, 1);
for i = 1:profileLen
    if any(isnan(data(i:i+subLen-1))) || any(isinf(data(i:i+subLen-1)))
        isSkip(i) = true;
    end
end
data(isnan(data)|isinf(data)) = 0;
    
%% preprocess for matrix profile
[dataFreq, data2Sum, dataSum, dataMean, data2Sig, dataSig] = ...
    massPre(data, dataLen, subLen);
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
    if isSkip(idx)
       continue 
    end
    query = data(idx:idx+subLen-1);
    if i == 1
        distanceProfile = mass(dataFreq, query, dataLen, subLen, ...
            data2Sum, dataSum, dataMean, data2Sig, dataSig);
        distanceProfile = abs(distanceProfile);
    else
        % replace with yan's method
        distanceProfile = mass(dataFreq, query, dataLen, subLen, ...
            data2Sum, dataSum, dataMean, data2Sig, dataSig);
        distanceProfile = abs(distanceProfile);
    end
    
    % apply skip zone
    distanceProfile(isSkip) = inf;
    
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
        
        % remove motif
        if exist('motifDataPlot', 'var')
            for j = 1:2
                delete(motifDataPlot(j));
            end
        end
        if exist('discordPlot', 'var')
            for j = 1:length(discordPlot)
                delete(discordPlot(j));
            end
        end
        if exist('motifMotifPlot', 'var')
            for j = 1:3
                for k = 1:2
                    for l = 1:length(motifMotifPlot{j, k})
                        delete(motifMotifPlot{j, k}(l));
                    end
                end
            end
        end
        
        % apply discard
        mainWindow = get(mainWindow.fig, 'userdata');
        discardIdx = mainWindow.discardIdx;
        matrixProfileTemp = matrixProfile;
        for j = 1:length(discardIdx)
            discardZoneStart = max(1, discardIdx(j)-exclusionZone);
            discardZoneEnd = min(profileLen, discardIdx(j)+exclusionZone);
            matrixProfileTemp(discardZoneStart:discardZoneEnd) = inf;
            matrixProfileTemp(abs(profileIndex - discardIdx(j)) < exclusionZone) = inf;
        end
        
        % compute matif
        motifIdxs = cell(3, 2);
        for j = 1:3
            [motifDistance, minIdx] = min(matrixProfileTemp);
            motifIdxs{j, 1} = sort([minIdx, profileIndex(minIdx)]);
            motifIdx = motifIdxs{j, 1}(1);
            motifQuery = data(motifIdx:motifIdx+subLen-1);
            
            % find neighbors
            motifDistanceProfile = mass(dataFreq, motifQuery, dataLen, subLen, ...
                data2Sum, dataSum, dataMean, data2Sig, dataSig);
            motifDistanceProfile = abs(motifDistanceProfile);
            motifDistanceProfile(motifDistanceProfile > motifDistance*radius) = inf;
            motifZoneStart = max(1, motifIdx-exclusionZone);
            motifZoneEnd = min(profileLen, motifIdx+exclusionZone);
            motifDistanceProfile(motifZoneStart:motifZoneEnd) = inf;
            motifIdx = motifIdxs{j, 1}(2);
            motifZoneStart = max(1, motifIdx-exclusionZone);
            motifZoneEnd = min(profileLen, motifIdx+exclusionZone);
            motifDistanceProfile(motifZoneStart:motifZoneEnd) = inf;
            motifDistanceProfile(isSkip) = inf;
            [distanceOrder, distanceIdxOrder] = sort(motifDistanceProfile, 'ascend');
            motifNeighbor = zeros(1, 10);
            for k = 1:10
                if isinf(distanceOrder(1)) || length(distanceOrder) < k
                    break;
                end
                motifNeighbor(k) = distanceIdxOrder(1);
                distanceOrder(1) = [];
                distanceIdxOrder(1) = [];
                distanceOrder(abs(distanceIdxOrder - motifNeighbor(k)) < exclusionZone) = [];
                distanceIdxOrder(abs(distanceIdxOrder - motifNeighbor(k)) < exclusionZone) = [];
            end
            motifNeighbor(motifNeighbor == 0) = [];
            motifIdxs{j, 2} = motifNeighbor;
            
            % remove found motif and their neighbor
            removeIdx = cell2mat(motifIdxs(j, :));
            for k = 1:length(removeIdx)
                removeZoneStart = max(1, removeIdx(k)-exclusionZone);
                removeZoneEnd = min(profileLen, removeIdx(k)+exclusionZone);
                matrixProfileTemp(removeZoneStart:removeZoneEnd) = inf;
            end
        end
        
        % plot motif on data
        motifColor = {'g', 'c'};
        motifDataPlot = zeros(2, 1);
        for j = 1:2
            motifPos = motifIdxs{1, 1}(j):motifIdxs{1, 1}(j)+subLen-1;
            hold(mainWindow.dataAx, 'on');
            motifDataPlot(j) = plot(motifPos, dataPlot(motifPos), ...
                motifColor{j}, 'parent', mainWindow.dataAx);
            hold(mainWindow.dataAx, 'off');
        end
        
        % plot motif's neighbor
        motifMotifPlot = cell(3, 2);
        neighborColor = 0.5*ones(1, 3);
        for j = 1:3
            motifMotifPlot{j, 2} = zeros(length(motifIdxs{j, 2}), 1);
            for k = 1:length(motifIdxs{j, 2})
                neighborPos = motifIdxs{j, 2}(k):motifIdxs{j, 2}(k)+subLen-1;
                if j == 1
                    hold(mainWindow.motif1Ax, 'on');
                    motifMotifPlot{j, 2}(k) = plot(1:subLen, zeroOneNorm(dataPlot(neighborPos)),...
                        'color', neighborColor, 'linewidth', 2, 'parent', mainWindow.motif1Ax);
                    hold(mainWindow.motif1Ax, 'off');
                elseif j == 2
                    hold(mainWindow.motif2Ax, 'on');
                    motifMotifPlot{j, 2}(k) = plot(1:subLen, zeroOneNorm(dataPlot(neighborPos)),...
                        'color', neighborColor, 'linewidth', 2, 'parent', mainWindow.motif2Ax);
                    hold(mainWindow.motif2Ax, 'off');
                elseif j == 3
                    hold(mainWindow.motif3Ax, 'on');
                    motifMotifPlot{j, 2}(k) = plot(1:subLen, zeroOneNorm(dataPlot(neighborPos)),...
                        'color', neighborColor, 'linewidth', 2, 'parent', mainWindow.motif3Ax);
                    hold(mainWindow.motif3Ax, 'off');
                end
            end
        end
        
        % plot motif on motif axis
        for j = 1:3
            motifMotifPlot{j, 1} = zeros(2, 1);
            for k = 1:2
                motifPos = motifIdxs{j, 1}(k):motifIdxs{j, 1}(k)+subLen-1;
                if j == 1
                    hold(mainWindow.motif1Ax, 'on');
                    set(mainWindow.motif1Text, 'string', ...
                        sprintf('The best-so-far 1st motifs are located at %d (green) and %d (cyan)', ...
                        motifIdxs{j, 1}(1), motifIdxs{j, 1}(2)));
                    motifMotifPlot{j, 1}(k) = plot(1:subLen, zeroOneNorm(dataPlot(motifPos)),...
                        motifColor{k}, 'parent', mainWindow.motif1Ax);
                    hold(mainWindow.motif1Ax, 'off');
                elseif j == 2
                    hold(mainWindow.motif2Ax, 'on');
                    set(mainWindow.motif2Text, 'string', ...
                        sprintf('The best-so-far 2nd motifs are located at %d (green) and %d (cyan)', ...
                        motifIdxs{j, 1}(1), motifIdxs{j, 1}(2)));
                    motifMotifPlot{j, 1}(k) = plot(1:subLen, zeroOneNorm(dataPlot(motifPos)),...
                        motifColor{k}, 'parent', mainWindow.motif2Ax);
                    hold(mainWindow.motif2Ax, 'off');
                elseif j == 3
                    hold(mainWindow.motif3Ax, 'on');
                    set(mainWindow.motif3Text, 'string', ...
                        sprintf('The best-so-far 3rd motifs are located at %d (green) and %d (cyan)', ...
                        motifIdxs{j, 1}(1), motifIdxs{j, 1}(2)));
                    motifMotifPlot{j, 1}(k) = plot(1:subLen, zeroOneNorm(dataPlot(motifPos)),...
                        motifColor{k}, 'parent', mainWindow.motif3Ax);
                    hold(mainWindow.motif3Ax, 'off');
                end
            end
        end
        
        % find discord
        matrixProfileTemp(isinf(matrixProfileTemp)) = -inf;
        [~, profileIdxOrder] = ...
            sort(matrixProfileTemp, 'descend');
        discordIdx = zeros(3, 1);
        for j = 1:3
            if length(profileIdxOrder) < j
                break
            end
            discordIdx(j) = profileIdxOrder(1);
            profileIdxOrder(1) = [];
            profileIdxOrder(abs(profileIdxOrder - discordIdx(j)) < exclusionZone) = [];
        end
        discordIdx(discordIdx == 0) = nan;
        
        % plot discord
        discordPlot = zeros(sum(~isnan(discordIdx)), 1);
        discordColor = {'b', 'r', 'g'};
        for j = 1:3
            if isnan(discordIdx(j))
                break;
            end
            discordPos = discordIdx(j):discordIdx(j)+subLen-1;
            hold(mainWindow.discordAx, 'on');
            discordPlot(j) = plot(1:subLen, zeroOneNorm(dataPlot(discordPos)),...
                discordColor{j}, 'parent', mainWindow.discordAx);
            hold(mainWindow.discordAx, 'off');
        end
        
        % update process
        set(mainWindow.dataText, 'string', ...
            sprintf('We are %.1f%% done: The input time series: The best-so-far motifs are color coded (see bottom panel)', i*100/profileLen));
        set(mainWindow.discordText, 'string', ...
            sprintf('The top three discords %d(blue), %d(red), %d(green)', discordIdx(1), discordIdx(2), discordIdx(3)));
        
        % show the figure
        if firstUpdate
            set(mainWindow.fig, 'userdata', mainWindow);
            set(mainWindow.fig, 'visible', 'on');
            firstUpdate = false;
        end
        
        % check for stop
        mainWindow = get(mainWindow.fig, 'userdata');
        mainWindow.motifIdxs = motifIdxs;
        set(mainWindow.fig, 'userdata', mainWindow);
        if mainWindow.stopping
            set(mainWindow.fig, 'name', 'UCR Interactive Matrix Profile Calculation (Stopped)');
            return;
        end
        if i == profileLen
            set(mainWindow.fig, 'name', 'UCR Interactive Matrix Profile Calculation (Completed)');
            set(mainWindow.discard1Btn, 'enable', 'off');
            set(mainWindow.discard2Btn, 'enable', 'off');
            set(mainWindow.discard3Btn, 'enable', 'off');
            set(mainWindow.stopBtn, 'enable', 'off');
            return;
        end
        
        % pause for plot and restart timer
        set(mainWindow.discard1Btn, 'enable', 'on');
        set(mainWindow.discard2Btn, 'enable', 'on');
        set(mainWindow.discard3Btn, 'enable', 'on');
        pause(0.01);
        timer = tic();
    end
end


function pushDiscardBtn(src, ~, btnNum)
mainWindowFig = get(src, 'parent');
mainWindow = get(mainWindowFig, 'userdata');
mainWindow.discardIdx = [mainWindow.discardIdx, ...
    mainWindow.motifIdxs{btnNum, 1}];
set(mainWindow.discard1Btn, 'enable', 'off');
set(mainWindow.discard2Btn, 'enable', 'off');
set(mainWindow.discard3Btn, 'enable', 'off');
set(mainWindow.fig, 'userdata', mainWindow);


function pushStopBtn(src, ~)
mainWindowFig = get(src, 'parent');
mainWindow = get(mainWindowFig, 'userdata');
mainWindow.stopping = true;
set(mainWindow.discard1Btn, 'enable', 'off');
set(mainWindow.discard2Btn, 'enable', 'off');
set(mainWindow.discard3Btn, 'enable', 'off');
set(src, 'enable', 'off');
set(mainWindow.fig, 'userdata', mainWindow);


%% The following two functions are modified from the code provided in the following URL
%  http://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html
function [dataFreq, data2Sum, dataSum, dataMean, data2Sig, dataSig] = ...
    massPre(data, dataLen, subLen)
data(dataLen+1:2*dataLen) = 0;
dataFreq = fft(data);
cum_sumx = cumsum(data);
cum_sumx2 =  cumsum(data.^2);
data2Sum = cum_sumx2(subLen:dataLen)-[0;cum_sumx2(1:dataLen-subLen)];
dataSum = cum_sumx(subLen:dataLen)-[0;cum_sumx(1:dataLen-subLen)];
dataMean = dataSum./subLen;
data2Sig = (data2Sum./subLen)-(dataMean.^2);
dataSig = sqrt(data2Sig);


function distanceProfile = mass(dataFreq, query, dataLen, subLen, ...
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
x = x - min(x(~isinf(x) & ~isnan(x)));
x = x / max(x(~isinf(x) & ~isnan(x)));


function mainResize(src, ~)
mainWindow = get(src, 'userdata');
figPosition = get(mainWindow.fig, 'position');
axGap = 38;
axesHeight = round((figPosition(4)-axGap*5-60)/6);
set(mainWindow.dataAx, 'position', ...
    [30, 5*axesHeight+5*axGap+30, figPosition(3)-60, axesHeight]);
set(mainWindow.profileAx, 'position', ...
    [30, 4*axesHeight+4*axGap+30, figPosition(3)-60, axesHeight]);
set(mainWindow.motif1Ax, 'position', ...
    [30, 3*axesHeight+3*axGap+30, figPosition(3)-160, axesHeight]);
set(mainWindow.motif2Ax, 'position', ...
    [30, 2*axesHeight+2*axGap+30, figPosition(3)-160, axesHeight]);
set(mainWindow.motif3Ax, 'position', ...
    [30, 1*axesHeight+1*axGap+30, figPosition(3)-160, axesHeight]);
set(mainWindow.discordAx, 'position', ...
    [30, 30, figPosition(3)-160, axesHeight]);
set(mainWindow.discard1Btn, 'position', ...
    [figPosition(3)-120, 3*axesHeight+3*axGap+30, 90, 20]);
set(mainWindow.discard2Btn, 'position', ...
    [figPosition(3)-120, 2*axesHeight+2*axGap+30, 90, 20]);
set(mainWindow.discard3Btn, 'position', ...
    [figPosition(3)-120, 1*axesHeight+1*axGap+30, 90, 20]);
set(mainWindow.stopBtn, 'position', ...
    [figPosition(3)-120, 30, 90, 20]);
set(mainWindow.dataText, 'position', ...
    [30, 6*axesHeight+5*axGap+30, figPosition(3)-60, 18]);
set(mainWindow.profileText, 'position', ...
    [30, 5*axesHeight+4*axGap+30, figPosition(3)-60, 18]);
set(mainWindow.motif1Text, 'position', ...
    [30, 4*axesHeight+3*axGap+30, figPosition(3)-160, 18]);
set(mainWindow.motif2Text, 'position', ...
    [30, 3*axesHeight+2*axGap+30, figPosition(3)-160, 18]);
set(mainWindow.motif3Text, 'position', ...
    [30, 2*axesHeight+1*axGap+30, figPosition(3)-160, 18]);
set(mainWindow.discordText, 'position', ...
    [30, 1*axesHeight+30, figPosition(3)-160, 18]);