% File: IntegerBasedStructure(x,q,delta,d)
%
% Goal: find the data sites located in each of the q^d blocks
%
% Inputs: x:               nXd matrix representing a set of data
%             q:               number of blocks in one direction
%             delta:         radius of PU subdomains
%             d:               space dimension
%
% Outputs: X_block: multiarray containing the indices of the data points located in k-th block
%
function [X_block] = IntegerBasedStructure(x,q,delta,d)
n = size(x,1); X_block = cell(q^d,1); k = 1:d-1; % Initialize
for i = 1:n % Find the indices of the data points located in k-th block
    idx = ceil(x(i,:)/delta); idx(idx == 0) = 1; index = sum((idx(k)-1)*q.^(d-k)) + idx(end); 
    X_block{index} = [X_block{index}; i];
end