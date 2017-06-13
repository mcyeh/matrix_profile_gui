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
excZoneLen = round(subLen * 0.5);
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
mainWindow.dataAx = axes('parent', mainWindow.fig, ...
    'units', 'pixels', 'xlim', [1, dataLen], 'xtick', [1, dataLen], ...
    'ylim', [-0.05, 1.05], 'ytick', [], 'ycolor', [1, 1, 1]);
mainWindow.profileAx = axes('parent', mainWindow.fig, ...
    'units', 'pixels', 'xlim', [1, dataLen], 'xtick', [1, dataLen], ...
    'ylim', [0, 2*sqrt(subLen)]);
mainWindow.discordAx = axes('parent', mainWindow.fig, ...
    'units', 'pixels', 'xlim', [1, subLen], 'xtick', [1, subLen], ...
    'ylim', [-0.05, 1.05], 'ytick', [], 'ycolor', [1, 1, 1]);
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
mainWindow.discordText = uicontrol('parent', mainWindow.fig, ...
    'style', 'text', 'string', '', 'fontsize', 10, ...
    'backgroundcolor', backColor, 'horizontalalignment', 'left');
mainWindow.motifAx = zeros(3, 1);
for i = 1:3
    mainWindow.motifAx(i) = axes('parent', mainWindow.fig, ...
        'units', 'pixels', 'xlim', [1, subLen], 'xtick', [1, subLen], ...
        'ylim', [-0.05, 1.05], 'ytick', [], 'ycolor', [1, 1, 1]);
end
mainWindow.motifText = zeros(3, 1);
for i = 1:3
    mainWindow.motifText(i) = uicontrol('parent', mainWindow.fig, ...
        'style', 'text', 'string', '', 'fontsize', 10, ...
        'backgroundcolor', backColor, 'horizontalalignment', 'left');
end
mainWindow.discardBtn = zeros(3, 1);
for i = 1:3
    mainWindow.discardBtn(i) = uicontrol('parent',mainWindow.fig, ...
        'style', 'pushbutton', 'string', 'Discard', 'fontsize', 10, ...
        'callback', @(src, cbdata) pushDiscardBtn(src, cbdata, i));
end

%% plot data
dataPlot = zeroOneNorm(data);
hold(mainWindow.dataAx, 'on');
plot(1:dataLen, dataPlot, 'r', 'parent', mainWindow.dataAx);
hold(mainWindow.dataAx, 'off');

%% locate nan and inf
proLen = dataLen - subLen + 1;
isSkip = false(proLen, 1);
for i = 1:proLen
    if any(isnan(data(i:i + subLen - 1))) || ...
            any(isinf(data(i:i + subLen - 1)))
        isSkip(i) = true;
    end
end
data(isnan(data) | isinf(data)) = 0;

%% preprocess for matrix profile
[dataFreq, dataMu, dataSig] = massPre(data, dataLen, subLen);
idxOrder = randperm(proLen);
matrixProfile = inf(proLen, 1);
profileIndex = zeros(proLen, 1);

%% color and text settings
txtTemp = {
    'The best-so-far 1st motifs are located at %d (green) and %d (cyan)';
    'The best-so-far 2nd motifs are located at %d (green) and %d (cyan)';
    'The best-so-far 3rd motifs are located at %d (green) and %d (cyan)';
    };
motifColor = {'g', 'c'};
discordColor = {'b', 'r', 'g'};
neighborColor = 0.5 * ones(1, 3);

%% iteratively plot
mainWindow.stopping = false;
mainWindow.discardIdx = [];
set(mainWindow.fig, 'userdata', mainWindow);
firstUpdate = true;
timer = tic();
for i = 1:proLen
    idx = idxOrder(i);
    if isSkip(idx)
        continue
    end
    drawnow;
    
    % compute the distance profile
    query = data(idx:idx+subLen-1);
    if i == 1
        distProfile = mass(dataFreq, query, dataLen, subLen, ...
            dataMu, dataSig, dataMu(idx), dataSig(idx));
        distProfile = real(distProfile);
        distProfile = sqrt(distProfile);
    else
        distProfile = mass(dataFreq, query, dataLen, subLen, ...
            dataMu, dataSig, dataMu(idx), dataSig(idx));
        distProfile = real(distProfile);
        distProfile = sqrt(distProfile);
    end
    
    % apply exclusion and skip zone
    distProfile(isSkip) = inf;
    excZoneStart = max(1, idx - excZoneLen);
    excZoneEnd = min(proLen, idx + excZoneLen);
    distProfile(excZoneStart:excZoneEnd) = inf;
    
    % update matrix profile
    updatePos = distProfile < matrixProfile;
    profileIndex(updatePos) = idx;
    matrixProfile(updatePos) = distProfile(updatePos);
    [matrixProfile(idx), profileIndex(idx)] = min(distProfile);
    
    % check update condition
    if toc(timer) < 1 && i ~= proLen
        continue;
    end
    
    % plot matrix profile
    if exist('prefilePlot', 'var')
        delete(prefilePlot);
    end
    hold(mainWindow.profileAx, 'on');
    prefilePlot = plot(1:proLen, matrixProfile, ...
        'b', 'parent', mainWindow.profileAx);
    hold(mainWindow.profileAx, 'off');
    
    % remove motif
    if exist('motifMarkPlot', 'var')
        for j = 1:2
            delete(motifMarkPlot(j));
        end
    end
    if exist('discordPlot', 'var')
        for j = 1:length(discordPlot)
            delete(discordPlot(j));
        end
    end
    if exist('motifPlot', 'var')
        for j = 1:3
            for k = 1:2
                for l = 1:length(motifPlot{j, k})
                    delete(motifPlot{j, k}(l));
                end
            end
        end
    end
    
    % apply discard
    mainWindow = get(mainWindow.fig, 'userdata');
    discardIdx = mainWindow.discardIdx;
    matrixProfileCur = matrixProfile;
    for j = 1:length(discardIdx)
        discardZoneStart = max(1, discardIdx(j)-excZoneLen);
        discardZoneEnd = min(proLen, discardIdx(j)+excZoneLen);
        matrixProfileCur(discardZoneStart:discardZoneEnd) = inf;
        matrixProfileCur(abs(profileIndex - discardIdx(j)) < ...
            excZoneLen) = inf;
    end
    
    % compute matif
    motifIdxs = cell(3, 2);
    for j = 1:3
        [motifDistance, minIdx] = min(matrixProfileCur);
        motifIdxs{j, 1} = sort([minIdx, profileIndex(minIdx)]);
        motifIdx = motifIdxs{j, 1}(1);
        motifQuery = data(motifIdx:motifIdx+subLen-1);
        
        % find neighbors
        motifDistProfile = mass(dataFreq, motifQuery, ...
            dataLen, subLen, dataMu, dataSig, ...
            dataMu(motifIdx), dataSig(motifIdx));
        motifDistProfile = abs(motifDistProfile);
        motifDistProfile(motifDistProfile > motifDistance*radius) = inf;
        motifZoneStart = max(1, motifIdx-excZoneLen);
        motifZoneEnd = min(proLen, motifIdx+excZoneLen);
        motifDistProfile(motifZoneStart:motifZoneEnd) = inf;
        motifIdx = motifIdxs{j, 1}(2);
        motifZoneStart = max(1, motifIdx-excZoneLen);
        motifZoneEnd = min(proLen, motifIdx+excZoneLen);
        motifDistProfile(motifZoneStart:motifZoneEnd) = inf;
        motifDistProfile(isSkip) = inf;
        [distanceOrder, distanceIdxOrder] = sort(motifDistProfile, 'ascend');
        motifNeighbor = zeros(1, 10);
        for k = 1:10
            if isinf(distanceOrder(1)) || length(distanceOrder) < k
                break;
            end
            motifNeighbor(k) = distanceIdxOrder(1);
            distanceOrder(1) = [];
            distanceIdxOrder(1) = [];
            distanceOrder(abs(distanceIdxOrder - motifNeighbor(k)) < excZoneLen) = [];
            distanceIdxOrder(abs(distanceIdxOrder - motifNeighbor(k)) < excZoneLen) = [];
        end
        motifNeighbor(motifNeighbor == 0) = [];
        motifIdxs{j, 2} = motifNeighbor;
        
        % remove found motif and their neighbor
        removeIdx = cell2mat(motifIdxs(j, :));
        for k = 1:length(removeIdx)
            removeZoneStart = max(1, removeIdx(k)-excZoneLen);
            removeZoneEnd = min(proLen, removeIdx(k)+excZoneLen);
            matrixProfileCur(removeZoneStart:removeZoneEnd) = inf;
        end
    end
    
    % plot motif on data
    motifMarkPlot = zeros(2, 1);
    for j = 1:2
        motifPos = motifIdxs{1, 1}(j):motifIdxs{1, 1}(j)+subLen-1;
        hold(mainWindow.dataAx, 'on');
        motifMarkPlot(j) = plot(motifPos, dataPlot(motifPos), ...
            motifColor{j}, 'parent', mainWindow.dataAx);
        hold(mainWindow.dataAx, 'off');
    end
    
    % plot motif's neighbor
    motifPlot = cell(3, 2);
    for j = 1:3
        motifPlot{j, 2} = zeros(length(motifIdxs{j, 2}), 1);
        for k = 1:length(motifIdxs{j, 2})
            neighborPos = motifIdxs{j, 2}(k):motifIdxs{j, 2}(k)+subLen-1;
            
            hold(mainWindow.motifAx(j), 'on');
            motifPlot{j, 2}(k) = plot(1:subLen, ...
                zeroOneNorm(dataPlot(neighborPos)),...
                'color', neighborColor, 'linewidth', 2, ...
                'parent', mainWindow.motifAx(j));
            hold(mainWindow.motifAx(j), 'off');
        end
    end
    
    % plot motif on motif axis
    for j = 1:3
        motifPlot{j, 1} = zeros(2, 1);
        for k = 1:2
            motifPos = motifIdxs{j, 1}(k):motifIdxs{j, 1}(k)+subLen-1;
            
            hold(mainWindow.motifAx(j), 'on');
            set(mainWindow.motifText(j), 'string', ...
                sprintf(txtTemp{j}, ...
                motifIdxs{j, 1}(1), motifIdxs{j, 1}(2)));
            motifPlot{j, 1}(k) = plot(1:subLen, ...
                zeroOneNorm(dataPlot(motifPos)),...
                motifColor{k}, 'parent', mainWindow.motifAx(j));
            hold(mainWindow.motifAx(j), 'off');
        end
    end
    
    % find discord
    matrixProfileCur(isinf(matrixProfileCur)) = -inf;
    [~, profileIdxOrder] = ...
        sort(matrixProfileCur, 'descend');
    discordIdx = zeros(3, 1);
    for j = 1:3
        if length(profileIdxOrder) < j
            break
        end
        discordIdx(j) = profileIdxOrder(1);
        profileIdxOrder(1) = [];
        profileIdxOrder(abs(profileIdxOrder - discordIdx(j)) < ...
            excZoneLen) = [];
    end
    discordIdx(discordIdx == 0) = nan;
    
    % plot discord
    discordPlot = zeros(sum(~isnan(discordIdx)), 1);
    for j = 1:3
        if isnan(discordIdx(j))
            break;
        end
        discordPos = discordIdx(j):discordIdx(j)+subLen-1;
        hold(mainWindow.discordAx, 'on');
        discordPlot(j) = plot(1:subLen, ...
            zeroOneNorm(dataPlot(discordPos)),...
            discordColor{j}, 'parent', mainWindow.discordAx);
        hold(mainWindow.discordAx, 'off');
    end
    
    % update process text
    set(mainWindow.dataText, 'string', ...
        sprintf(['We are %.1f%% done: The input time series: ', ...
        'The best-so-far motifs are color coded (see bottom panel)'], ...
        i*100/proLen));
    set(mainWindow.discordText, 'string', ...
        sprintf(['The top three discords ', ...
        '%d(blue), %d(red), %d(green)'], ...
        discordIdx(1), discordIdx(2), discordIdx(3)));
    
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
    if i == proLen
        set(mainWindow.fig, 'name', ...
            'UCR Interactive Matrix Profile Calculation (Completed)');
    elseif mainWindow.stopping
        set(mainWindow.fig, 'name', ...
            'UCR Interactive Matrix Profile Calculation (Stopped)');
    end
    if i == proLen || mainWindow.stopping
        for j = 1:3
            set(mainWindow.discardBtn(j), 'enable', 'off');
        end
        set(mainWindow.stopBtn, 'enable', 'off');
        return;
    end
    
    % restart timer
    for j = 1:3
        set(mainWindow.discardBtn(j), 'enable', 'on');
    end
    timer = tic();
end


function pushDiscardBtn(src, ~, btnNum)
mainWindowFig = get(src, 'parent');
mainWindow = get(mainWindowFig, 'userdata');
mainWindow.discardIdx = [mainWindow.discardIdx, ...
    mainWindow.motifIdxs{btnNum, 1}];
for i = 1:3
    set(mainWindow.discardBtn(i), 'enable', 'off');
end
set(mainWindow.fig, 'userdata', mainWindow);


function pushStopBtn(src, ~)
mainWindowFig = get(src, 'parent');
mainWindow = get(mainWindowFig, 'userdata');
mainWindow.stopping = true;
for i = 1:3
    set(mainWindow.discardBtn(i), 'enable', 'off');
end
set(src, 'enable', 'off');
set(mainWindow.fig, 'userdata', mainWindow);


%% The following two functions are modified from the code provided in the following URL
%  http://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html
function [dataFreq, dataMu, dataSig] = massPre(data, dataLen, subLen)
data(dataLen + 1:(subLen + dataLen)) = 0;
dataFreq = fft(data);
dataCumsum = cumsum(data);
data2Cumsum =  cumsum(data .^ 2);
data2Sum = data2Cumsum(subLen:dataLen) - ...
    [0; data2Cumsum(1:dataLen - subLen)];
dataSum = dataCumsum(subLen:dataLen) - ...
    [0; dataCumsum(1:dataLen - subLen)];
dataMu = dataSum ./ subLen;
data2Sig = (data2Sum ./ subLen) - (dataMu .^ 2);
dataSig = sqrt(data2Sig);


function distProfile = mass(dataFreq, query, ...
    dataLen, subLen, dataMu, dataSig, queryMu, querySig)
query = query(end:-1:1);
query(subLen+1:(subLen+dataLen)) = 0;
queryFreq = fft(query);
productFreq = dataFreq .* queryFreq;
product = ifft(productFreq);
distProfile = 2 * (subLen - ...
    (product(subLen:dataLen) - subLen * dataMu * queryMu) ./ ...
    (dataSig * querySig));


function x = zeroOneNorm(x)
x = x - min(x(~isinf(x) & ~isnan(x)));
x = x / max(x(~isinf(x) & ~isnan(x)));


function mainResize(src, ~)
mainWindow = get(src, 'userdata');
figPosition = get(mainWindow.fig, 'position');
axGap = 38;
axesHeight = round((figPosition(4) - axGap * 5 - 60) / 6);
set(mainWindow.dataAx, 'position', ...
    [30, 5 * axesHeight+5 * axGap + 30, figPosition(3) - 60, axesHeight]);
set(mainWindow.profileAx, 'position', ...
    [30, 4 * axesHeight+4 * axGap + 30, figPosition(3) - 60, axesHeight]);
set(mainWindow.discordAx, 'position', ...
    [30, 30, figPosition(3) - 160, axesHeight]);
set(mainWindow.stopBtn, 'position', ...
    [figPosition(3) - 120, 30, 90, 20]);
set(mainWindow.dataText, 'position', ...
    [30, 6 * axesHeight + 5 * axGap + 30, figPosition(3) - 60, 18]);
set(mainWindow.profileText, 'position', ...
    [30, 5 * axesHeight + 4 * axGap + 30, figPosition(3) - 60, 18]);
set(mainWindow.discordText, 'position', ...
    [30, 1 * axesHeight + 30, figPosition(3) - 160, 18]);
for i = 1:3
    set(mainWindow.motifAx(i), 'position', ...
        [30, (4 - i) * axesHeight + (4 - i) * axGap + 30, ...
        figPosition(3) - 160, axesHeight]);
end
for i = 1:3
    set(mainWindow.motifText(i), 'position', ...
        [30, (5 - i) * axesHeight + (4 - i) * axGap + 30, ...
        figPosition(3) - 160, 18]);
end
for i = 1:3
    set(mainWindow.discardBtn(i), 'position', ...
        [figPosition(3) - 120, ...
        (4 - i) * axesHeight + (4 - i) * axGap + 30, 90, 20]);
end