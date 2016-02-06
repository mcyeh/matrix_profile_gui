clear;
clc;

%%
subLen = 2^8;
dataLen = 2^14;

%%
dataLen = dataLen-subLen*6;
randonLen = round(dataLen/7);

%%
pattern1 = sin(0:(2*pi)/subLen:2*pi)+2*tanh(0:(2*pi)/subLen:2*pi);
pattern1 = pattern1(1:subLen);
pattern1 = pattern1 - min(pattern1);
pattern1 = pattern1 / max(pattern1);
pattern1 = pattern1 * 2 - 1;

pattern2 = sin(2*pi:-(2*pi)/subLen:0)+4*tanh(2*pi:-(2*pi)/subLen:0);
pattern2 = pattern2(1:subLen);
pattern2 = pattern2 - min(pattern2);
pattern2 = pattern2 / max(pattern2);
pattern2 = pattern2 * 2 - 1;

pattern2Anti = pattern2;
pattern2Anti = pattern2Anti - min(pattern2Anti);
pattern2Anti = pattern2Anti / max(pattern2Anti);
pattern2Anti = max(pattern2Anti) - pattern2Anti;
pattern2Anti = pattern2Anti * 2 - 1;

%%
data = cell(1, 13);
data{2} = pattern1+randn(size(pattern1))*0.01;
data{4} = pattern1+randn(size(pattern1))*0.01;
data{6} = pattern1+randn(size(pattern1))*0.01;

% data{2} = pattern2+randn(size(pattern2))*0.01;
% data{4} = pattern2+randn(size(pattern2))*0.01;
% data{6} = pattern2+randn(size(pattern2))*0.01;

data{8} = pattern2+randn(size(pattern2))*0.01;
data{10} = pattern2+randn(size(pattern2))*0.01;

data{12} = pattern2Anti+randn(size(pattern2Anti))*0.01;

%%
% data{1} = cumsum(randn(1, randonLen)*0.1);
% data{3} = cumsum(randn(1, randonLen)*0.1);
% data{5} = cumsum(randn(1, randonLen)*0.1);
% data{7} = cumsum(randn(1, randonLen)*0.1);
% data{9} = cumsum(randn(1, randonLen)*0.1);
% data{11} = cumsum(randn(1, dataLen-5*randonLen)*0.1);
data{1} = randn(1, randonLen)*0.1;
data{3} = randn(1, randonLen)*0.1;
data{5} = randn(1, randonLen)*0.1;
data{7} = randn(1, randonLen)*0.1;
data{9} = randn(1, randonLen)*0.1;
data{11} = randn(1, randonLen)*0.1;
data{13} = randn(1, dataLen-6*randonLen)*0.1;
patPos = [randonLen + 1, 2*randonLen + subLen + 1, 3*randonLen + 2*subLen + 1, 4*randonLen + 3*subLen + 1, 5*randonLen + 4*subLen + 1, 6*randonLen + 5*subLen + 1];

%%
data = cell2mat(data);

%%
figure;
hold on;
plot(1:length(data), data, 'b');
for i = 1:length(patPos)
    if i <= 3
        plot(patPos(i):patPos(i)+subLen-1, ...
            data(patPos(i):patPos(i)+subLen-1), 'r');
    elseif i <= 5
        plot(patPos(i):patPos(i)+subLen-1, ...
            data(patPos(i):patPos(i)+subLen-1), 'g');
    else
        plot(patPos(i):patPos(i)+subLen-1, ...
            data(patPos(i):patPos(i)+subLen-1), 'c');
    end
end
hold off;
xlim(gca, [1, dataLen]);

%%
save('testData', 'data', 'subLen')