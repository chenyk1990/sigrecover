function [ D1 ] = sr_mf(D,nfw,ifb,axis)
%sr_MF: median filter along first or second axis for 2D profile
%  IN   D:   	intput data 
%       nfw:    window size 
%       ifb:    if use padded boundary (if not, zero will be padded)
%       axis:    temporal sampling interval
%      
%  OUT   D1:  	output data
% 
%  Copyright (C) 2014 The University of Texas at Austin
%  Copyright (C) 2014 Yangkang Chen
%
% Example: dsp/test_sr_meanf_repeat.m
%
%
% References
% Huang et al., 2021, Erratic noise suppression using iterative structure-oriented space-varying median filtering with sparsity constraint, Geophysical Prospecting, 69, 101-121.
% Chen et al., 2020, Deblending of simultaneous-source data using a structure-oriented space-varying median filter, Geophysical Journal International, 222, 1805–1823.
% Gan et al., 2016, Separation of simultaneous sources using a structural-oriented median filter in the flattened dimension, Computers & Geosciences, 86, 46-54.
% Chen, Y., 2015, Deblending using a space-varying median filter, Exploration Geophysics, 46, 332-341.

if nargin==0
 error('Input data must be provided!');
end

if nargin==1
 nfw=7;
 ifb=1;
 axis=2;
end;

if nargin==2
 ifb=1;
 axis=2;    
end

% nfw should be odd
if mod(nfw,2)==0
    nfw=nfw+1;
end

if axis==2
   D=D.'; 
end

[n1,n2]=size(D);
nfw2=(nfw-1)/2;

if ifb==1
    D=[flipud(D(1:nfw2,:));D;flipud(D(n1-nfw2+1:n1,:))];
else
    D=[zeros(nfw2,n2);D;zeros(nfw2,n2)];    
end

% output data
D1=zeros(n1,n2);
for i2=1:n2
   for i1=1:n1
      D1(i1,i2)=median(D(i1:i1+nfw-1,i2)); 
   end 
end
if axis==2
    D1=D1.';
end
return
