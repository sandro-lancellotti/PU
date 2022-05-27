% File: IntegerBasedContainingQuery(tilde_x,q,delta,d)
%
% Goal: script that given a subdomain centre returns the index of the square block containing the 
%          subdomain centre
%
% Inputs: tilde_x:         subdomain centre
%             q:                 number of blocks in one direction
%             delta:           radius of the PU subdomains
%             d:                 space dimension
%
% Outputs: k: the index of the block containing the subdomain centre
%
function [k] = IntegerBasedContainingQuery(tilde_x,q,delta,d)
k_l = ceil(tilde_x/delta); l = 1:d-1; k_l(k_l == 0) = 1; k = sum((k_l(l)-1)*q.^(d-l)) + k_l(end);