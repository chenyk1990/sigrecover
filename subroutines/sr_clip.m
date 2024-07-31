function [dout]=sr_clip(din,mmin,mmax)
% sr_clip: clip data
% Yangkang Chen
% Aug, 2, 2020
% 
% EXAMPLE: 
% d=magic(8)
% sr_clip(d,5,40)
% 
[n1,n2,n3,n4,n5]=size(din);

dout=din;
tmp=find(dout>mmax);
dout(tmp)=mmax*ones(size(tmp));

tmp=find(dout<mmin);
dout(tmp)=mmin*ones(size(tmp));

dout=reshape(dout,n1,n2,n3,n4,n5);
return