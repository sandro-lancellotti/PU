% -------------------------------------------------------------------------------- %
%      This software is licensed by Creative Commons BF-NC-SA:    %
%  http://creativecommons.org/licenses/by-nc-sa/3.0/it/legalcode %
% -------------------------------------------------------------------------------- %
%
% File: PU(d,x,bar_x,m_d,phi,w,f,epsilon)
%
% Goal: script that performs partition of unity
%
% INPUTS:  d:              space dimension
%               x:               nXd matrix representing a set of n interpolation data
%               bar_x:        sXd matrix representing a set of s evaluation data
%               m_d:          number of PU subdomains in one direction
%               phi:            radial basis function
%               w:             weight function
%               f:               the function values
%               epsilon:      the shape parameter
%
% Outputs:  Pf:           sXd matrix representing the PU fit
%
% Calls on: IntegerBasedStructure,
%                IntegerBasedNeighbourhood,
%                IntegerBasedRangeSearch,
%                IntegerBasedContainingQuery,
%                DistanceMatrix, MakeSDGrid
%
% Remark:   DistanceMatrix and MakeSDGrid come from the book:
%           [G. bar_x. Fasshauer, Meshfree approximation methods with
%           Matlab, World Scientific, Singapore, 2007].
%                 IntegerBasedStructure,
%                 IntegerBasedNeighbourhood,
%                 IntegerBasedRangeSearch,
%                 IntegerBasedContainingQuery, are  by
%           R. Cavoretto, A. De Rossi, bar_x. Perracchione, are avilable at
%           available at: http://hdl.handle.net/2318/1559094
%
% -------------------------------------------------------------------------------- %
function [Pf] = PU(d,x,bar_x,m_d,phi,w,f,epsilon)
tilde_x = MakeSDGrid(d,m_d); % Create m_d^d equally spaced PU centres
delta = sqrt(2)/m_d;  supp = 1/delta; % Define the PU radius and the parameter for the weight functions
m = size(tilde_x,1); s = size(bar_x,1); Pf = zeros(s,1);  % Initialize and compute the Shepard matrix
DM_eval = DistanceMatrix(bar_x,tilde_x); W = w(supp,DM_eval); W = spdiags(1./(W*ones(m,1)),0,s,s)*W;
q = ceil(1./delta); % Parameter for the integer-based partitioning structure
% Build the partitioning structure for interpolation and evaluation data
X_block = IntegerBasedStructure(x,q,delta,d); bar_X_block = IntegerBasedStructure(bar_x,q,delta,d);
for j = 1:m % Loop over subdomains
    k = IntegerBasedContainingQuery(tilde_x(j,:),q,delta,d); % Find the box with the j-th PU centre
    % Find the interpolation data located in the j-th subdomain
    [X_NeigBlock, idx_X_NeigBlock] = IntegerBasedNeighbourhood(x,X_block,k,q,d);
    n_j = IntegerBasedRangeSearch(tilde_x(j,:),delta,X_NeigBlock,idx_X_NeigBlock);
    % Find the evaluation data located in the j-th subdomain
    [bar_X_NeigBlock, idx_bar_X_NeigBlock] = IntegerBasedNeighbourhood(bar_x,bar_X_block,k,q,d);
    s_j = IntegerBasedRangeSearch(tilde_x(j,:),delta,bar_X_NeigBlock,idx_bar_X_NeigBlock);
    if  (~isempty(s_j)) &&  (~isempty(n_j))
        DM_data = DistanceMatrix(x(n_j,:),x(n_j,:)); K_j = phi(epsilon,DM_data); % Interpolation matrix
        c_j = K_j\f(n_j);  % Compute the interpolation coefficients
        DM_eval = DistanceMatrix(bar_x(s_j,:),x(n_j,:)); % Compute the evaluation matrix
        K_eval = phi (epsilon, DM_eval); P_fj = K_eval*c_j; % Compute the the local fit
        Pf(s_j) = Pf(s_j) + P_fj.*W(s_j,j); % Accumulate the global fit
    end
end