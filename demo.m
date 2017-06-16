clear;
clc;

%% load test data
load('testData.mat');

matrixProfileTime = tic();
[matrixProfile, profileIndex, motifIndex, discordIndex] = ...
    interactiveMatrixProfile(data, subLen);
matrixProfileTime = toc(matrixProfileTime);

%% annotation vector framework demo
load('test_av_data');
av = make_AV_complexity(data, subLen);
tic;[matrixProfile, profileIndex, motifIndex, discordIndex] = ...
    interactiveMatrixProfile(data, subLen, av); toc; 