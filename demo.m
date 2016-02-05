clear;
clc;

%% load test data
load('testData.mat');


[matrixProfile, profileIndex, motifIndex] = interactiveMatrixProfile(data, subLen);