function [ X ] = sr_patch( A,mode,l1,l2,s1,s2 )
%decompose the image into patches:
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

[n1,n2]=size(A);

if mode==1 %possible for other patching options

    if nargin==2
        l1=8;l2=8;s1=4;s2=4;
    end

    if nargin==4
        s1=round(l1/2);s2=round(l2/2);
    end

    tmp=mod(n1-l1,s1);
    if tmp~=0
        A=[A;zeros(s1-tmp,n2)];
    end
    tmp=mod(n2-l2,s2);
    if tmp~=0
        A=[A,zeros(size(A,1),s2-tmp)];
    end



    [N1,N2]=size(A);
    X=[];
    for i1=1:s1:N1-l1+1
        for i2=1:s2:N2-l2+1
            tmp=reshape(A(i1:i1+l1-1,i2:i2+l2-1),l1*l2,1);
            X=[X,tmp];
        end
    end

end




end

