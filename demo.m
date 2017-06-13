clear;
clc;

%% load test data
load('testData.mat');

matrixProfileTime = tic();
[matrixProfile, profileIndex, motifIndex, discordIndex] = ...
    interactiveMatrixProfile(data, subLen);
matrixProfileTime = toc(matrixProfileTime);