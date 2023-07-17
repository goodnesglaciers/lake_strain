function [NewYr,NewX,NewCov]=FormMatrix(yr,x,cov,osite)

% FormMatrix
% [NewYr,NewX,NewCov]=FormMatrix(yr,x,cov,osite)
%
% Tools for shurinking matrix according to the number of data
% 
% June 11, 1997, Yosuke Aoki

NewYr=yr(:,1:osite);
NewX=x(:,1:osite);
NewCov=cov(1:3*osite,1:3*osite);