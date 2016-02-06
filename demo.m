clear;
clc;

%% load test data
load('testData.mat');


[matrixProfile, profileIndex, motifIndex, discordIndex] = interactiveMatrixProfile(data, subLen);