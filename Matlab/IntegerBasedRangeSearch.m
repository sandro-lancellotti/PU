% File: IntegerBasedRangeSearch(tilde_x,delta,X_NeigBlock,idx_X_NeigBlock)
%
% Goal: find the data sites located in a given subdomain and the distances between the 
%          subdomain centre and data sites
%
% Inputs:  tilde_x:                 subdomain centre
%             delta:                    radius of PU subdomains
%             X_NeigBlock:        nXd matrix representing a set of n points
%             idx_X_NeigBlock:  vector containing the indices of the data points located in the k-th block 
%                                         and in the neighbouring blocks
%
% Outputs: n_j:  vector containing the indices of the data points located in a given PU subdomain
%
function [n_j] = IntegerBasedRangeSearch(tilde_x,delta,X_NeigBlock,idx_X_NeigBlock)
N = size(X_NeigBlock,1); n_j = []; % Initialize
for i = 1:N % Compute distances between the data sites and the centre
    dist1(i) = norm(tilde_x - X_NeigBlock(i,:));
end
if N > 0 % Use a sort procedure to order distances
    [sort_dist,IX] = sort(dist1); N1 = size(sort_dist,2); j1 = 1; j2 = 1; %Initialize
    while (j2 <= N1) && (sort_dist(j2) <= delta) % Find the data sites located in the given subdomain
        n_j(j1) = idx_X_NeigBlock(IX(j2)); j1 = j1 + 1; j2 = j2 + 1;
    end
end