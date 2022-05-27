% File: IntegerBasedNeighbourhood(X,X_block,k,q,d)
%
% Goal: script that finds the neighbouring blocks
%
% Inputs:  x:           nXd matrix representing a set of n data sites
%              X_block: the integer-based data structure
%              k:           the k-th block cointaining the subdomain centre
%              q:           number of blocks in one direction
%              d:           space dimension
%
% Outputs:  X_NeigBlock:         points lying in the k-th neighbourhood
%                 idx_X_NeigBlock:   indices of points lying in the k-th neighbourhood
%
%-------------------------------------------------------------------------%
function [X_NeigBlock, idx_X_NeigBlock] = IntegerBasedNeighbourhood(x,X_block,k,q,d)
neigh = []; l = d-1; index = k; % Initialize
while l  > 0 % Find neighbouring blocks
     neigh = [k+q.^l,k-q.^l];
    if l - 1 > 0
        neigh = [neigh,neigh+q.^(l-1),neigh-q.^(l-1)];
    end
    l = l - 1;
end
k2 = 1; neighplus = []; neighminus = []; % Initialize
for i = 1:length(neigh)
    neighplus(k2) = neigh(i) + 1; neighminus(k2) = neigh(i) - 1; k2 = k2 + 1;
end
neigh = [neigh,k+1,k-1,neighplus,neighminus]; % Reduce the number of blocks for border blocks
j = find(neigh > 0 & neigh <= q^d); index = [index; neigh(j)']; X_NeigBlock = []; idx_X_NeigBlock = []; 
for p = 1:length(index)
    X_NeigBlock = [X_NeigBlock;x(X_block{index(p)},:)]; 
    idx_X_NeigBlock = [idx_X_NeigBlock;X_block{index(p)}];
end