% DeJong_f4_noise.m
% De Jong's f4 function, ND, includes a normally distributed noise term
%
% f(x) = sum( [1:N].*(in.^4), 2) + randn(tlen,1)
%
% x = N element row vector containing [ x0, x1,..., xN ]
%   each row is processed independently,
%   you can feed in matrices of timeXN no prob
%
% example: cost = DeJong_f4([1,2;3,4;5,6])
% note minimum @ x= all zeros

% Brian Birge
% Rev 1.0
% 9/12/04
function [out]=DeJong_f4_noise(in)
 persistent D tlen d randterm

% this speeds routine up a lot, if called from PSO these won't change from
% call to call (repmat is a cpu hog)
 Dx=length(in(1,:));
 tlenx=length(in(:,1));
 if isempty(D) | D~=Dx | tlen~=tlenx
   D=Dx; % dimension of prob
   tlen=tlenx; % how many separate states
   d=repmat([1:D],tlen,1); % needed to vectorize this
  % makes sure the function stays the same for the whole PSO run
   randterm=randn(tlen,1); 
 end
 out = sum( d.*(in.^4), 2) + randterm;