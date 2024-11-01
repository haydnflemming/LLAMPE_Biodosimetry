% -- DATA MIXER -----------------------------------------------------------
% Written by Cristian Ciobanu
%
% Given a 'set' that is a cell array consisting of data matrices (e.g.,
% 'set' = {control, dosed1, dosed2, dosed5}), this function will shuffle
% each matrix, then combine 70% of each matrix into a single 'train'
% matrix, and the remaining 30% of each matrix into a single 'test' matrix.

function [train, test] = combine(set)
    global C;
    C = {};
    train = [];
    test = [];
    numRows = 0;
    
    for i = set
        numRows = numRows + size(i{1}, 1);
    end
    
    for i = 1:max(size(set))
        [C{i, 1}, C{i, 2}] = shuffleSet(set{i});
        train = [train; C{i, 1}];
        test = [test; C{i, 2}];
    end
end

function [shuffledTrain, shuffledTest] = shuffleSet(X)
    numTrain = round(size(X, 1) * 0.7);
    shuffleX = X(randperm(size(X, 1)), :);
    shuffledTrain = shuffleX(1:numTrain, :);
    shuffledTest = shuffleX(numTrain + 1:end, :);
end