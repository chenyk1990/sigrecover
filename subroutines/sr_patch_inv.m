function [ A ] = sr_patch_inv( X,mode,n1,n2,l1,l2,s1,s2 )
% insert patches into the image
%  
% by Yangkang Chen
% Oct, 2017
%
% Input: 
%   D: input image
%   mode: patching mode
%   l1: first patch size
%   l2: second patch size
%   s1: first shifting size
%   s2: second shifting size
%   
% Output:
%   X: patches
% 
% Modified on Dec 12, 2018 (the edge issue, arbitrary size for the matrix)
% 		      Dec 31, 2018 (tmp1=mod(n1,l1) -> tmp1=mod(n1-l1,s1))
% 

%% patch size l1*l2
%l1=8;l2=8;
%

if mode==1 %possible for other patching options

if nargin==4
   l1=8;l2=8;s1=4;s2=4; 
end

if nargin==6
   s1=round(l1/2);s2=round(l2/2);
end

tmp1=mod(n1-l1,s1);
tmp2=mod(n2-l2,s2);
if tmp1~=0 && tmp2~=0
   A=zeros(n1+s1-tmp1,n2+s2-tmp2); 
   mask=zeros(n1+s1-tmp1,n2+s2-tmp2); 
end

if tmp1~=0 && tmp2==0
   A=zeros(n1+s1-tmp1,n2); 
   mask=zeros(n1+s1-tmp1,n2);
end

if tmp1==0 && tmp2~=0
   A=zeros(n1,n2+s2-tmp2);   
   mask=zeros(n1,n2+s2-tmp2);   
end

if tmp1==0 && tmp2==0
   A=zeros(n1,n2); 
   mask=zeros(n1,n2);
end

[N1,N2]=size(A);
id=0;
for i1=1:s1:N1-l1+1
    for i2=1:s2:N2-l2+1
        id=id+1;
%         [i1,i2,id]
        A(i1:i1+l1-1,i2:i2+l2-1)=A(i1:i1+l1-1,i2:i2+l2-1)+reshape(X(:,id),l1,l2);
        mask(i1:i1+l1-1,i2:i2+l2-1)=mask(i1:i1+l1-1,i2:i2+l2-1)+ones(l1,l2);
    end
end
A=A./mask; 
 
A=A(1:n1,1:n2);
end







end

