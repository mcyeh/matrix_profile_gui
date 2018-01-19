clear;
clc;

%% load test data
load('testData.mat');

%%
% without constrain, the seoncd and third pattern are the motif pair
matrixProfileTime1 = tic();
[matrixProfile1, profileIndex1, motifIndex1, discordIndex1] = ...
    interactiveMatrixProfile(data, subLen, inf); 
matrixProfileTime1 = toc(matrixProfileTime1);

%%
% by forcing the algorithm to cross a certain point, the first and third 
% pattern become the motif pair
matrixProfileTime2 = tic();
[matrixProfile2, profileIndex2, motifIndex2, discordIndex2] = ...
    interactiveMatrixProfile(data, subLen, 2800);
matrixProfileTime2 = toc(matrixProfileTime2);

fprintf('For the first part, the index of neighbor range from %d to %d\n', ...
    min(unique(profileIndex2(1:2800 - subLen + 1))), ...
    max(unique(profileIndex2(1:2800 - subLen + 1))));

fprintf('For the second part, the index of neighbor range from %d to %d\n', ...
    min(unique(profileIndex2(2800 + 1:end))), ...
    max(unique(profileIndex2(2800 + 1:end))));