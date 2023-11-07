clc;clear;close all;
addpath(genpath('../matseisdl'));

% DEMO script for 3D denoising and reconstruction based on SGK
% By Yangkang Chen
% March, 2020
%
% Key reference:
% Chen, Y., 2017, Fast dictionary learning for noise attenuation of multidimensional seismic data, Geophysical Journal International, 209, 21-31.
% Chen, Y., W. Huang, D. Zhang, and W. Chen, 2016, An open-source Matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
% 
% More References:
% Chen, Y., S. Fomel, 2015, Random noise attenuation using local signal-and-noise orthogonalization, Geophysics, 80, WD1-WD9.
% Chen, Y., J. Ma, and S. Fomel, 2016, Double-sparsity dictionary for seismic noise attenuation, Geophysics, 81, V17-V30.
% Siahsar, M. A. N., Gholtashi, S., Kahoo, A. R., W. Chen, and Y. Chen, 2017, Data-driven multi-task sparse dictionary learning for noise attenuation of 3D seismic data, Geophysics, 82, V385-V396.
% Siahsar, M. A. N., V. Abolghasemi, and Y. Chen, 2017, Simultaneous denoising and interpolation of 2D seismic data using data-driven non-negative dictionary learning, Signal Processing, 141, 309-321.
% Chen, Y., M. Zhang, M. Bai, and W. Chen, 2019, Improving the signal-to-noise ratio of seismological datasets by unsupervised machine learning, Seismological Research Letters, 90, 1552-1564.
% Chen, Y., S. Zu, W. Chen, M. Zhang, and Z. Guan, 2019, Learning the blending spikes using sparse dictionaries, Geophysical Journal International, 218, 1379?1397. 
% Wang, H., Q. Zhang, G. Zhang, J. Fang, and Y. Chen, 2020, Self-training and learning the waveform features of microseismic data using an adaptive dictionary, Geophysics, 85, KS51?KS61.
% Zu, S., H. Zhou, R. Wu, M. Jiang, and Y. Chen, 2019, Dictionary learning based on dip patch selection training for random noise attenuation, Geophysics, 84, V169?V183.
% Zu, S., H. Zhou, R. Wu, and Y. Chen, 2019, Hybrid-sparsity constrained dictionary learning for iterative deblending of extremely noisy simultaneous-source data, IEEE Transactions on Geoscience and Remote Sensing, 57, 2249-2262.
% etc. 

%% create data
a1=zeros(300,20);
[n,m]=size(a1);
a3=a1;
a4=a1;

k=0;
a=0.1;
b=1;
for t=-0.055:0.002:0.055
    k=k+1;
    b1(k)=(1-2*(pi*30*t).^2).*exp(-(pi*30*t).^2);
    b2(k)=(1-2*(pi*40*t).^2).*exp(-(pi*40*t).^2);
    b3(k)=(1-2*(pi*40*t).^2).*exp(-(pi*40*t).^2);
    b4(k)=(1-2*(pi*30*t).^2).*exp(-(pi*30*t).^2);
end
for i=1:m
  t1(i)=round(140);
  t3(i)=round(-6*i+180);
  t4(i)=round(6*i+10);
  a1(t1(i):t1(i)+k-1,i)=b1; 
  a3(t3(i):t3(i)+k-1,i)=b1; 
  a4(t4(i):t4(i)+k-1,i)=b1;
end

temp=a1(1:300,:)+a3(1:300,:)+a4(1:300,:);
for j=1:20
    a4=zeros(300,20);
    for i=1:m
  t4(i)=round(6*i+10+3*j); 
  a4(t4(i):t4(i)+k-1,i)=b1;
  
  t1(i)=round(140-2*j);
  a1(t1(i):t1(i)+k-1,i)=b1;
    end
    shot(:,:,j)=a1(1:300,:)+a3(1:300,:)+a4(1:300,:);
end
plane3d=shot;
dc=dl_scale(plane3d,3);
dc=dc(51:225,:,:);

%% select 1 slice from the 3D data as a 2D example
dc=dc(:,:,10);

%% adding noise
randn('state',201314);
var=0.1;
dn=dc+var*randn(size(dc));

% decimate
[nt,nx,ny]=size(dc);
ratio=0.65;
mask=dl_genmask(reshape(dc,nt,nx*ny),ratio,'c',201415);
mask=reshape(mask,nt,nx,ny);
d0=dn.*mask;
figure;imagesc([dc,dn,d0]);%colormap(seis);

%% simultaneous denoising and reconstruction
% Linear initialization (not useful ?)
% d1=dl_dlrecon_init(d0,mask,'nearest');
% % d1=InpaintingInterp2(d0,mask,'nearest');
% figure;imagesc([dc,d0,d1]);colormap(seis);


%% SGK
param=struct('T',3,'niter',10,'mode',1,'K',64);
mode=1;l1=6;l2=4;l3=1;s1=2;s2=2;s3=1;perc=2;Niter=12; %
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing for noisy data
% a=ones(Niter,1);
d1=dl_sgk_recon(d0,mask,mode,[l1,l2,l3],[s1,s2,s3],perc,Niter,a,param);
figure;imagesc([dc,dn,d0,d1]);%colormap(seis);
dl_snr(dc,d1)

%% KSVD
% param=struct('T',2,'niter',10,'mode',1,'K',64);
% mode=1;l1=4;l2=4;l3=4;s1=2;s2=2;s3=1;perc=2;Niter=10; %
% a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing for noisy data
% a=ones(Niter,1);
d2=dl_ksvd_recon(d0,mask,mode,[l1,l2,l3],[s1,s2,s3],perc,Niter,a,param);
figure;imagesc([dc,dn,d0,d2]);%colormap(seis);
dl_snr(dc,d2)

dl_snr(dc,dn) %
dl_snr(dc,d0) %
dl_snr(dc,d1) %
dl_snr(dc,d2) %


%% benchmark with DRR (https://github.com/chenyk1990/MATdrr)
addpath(genpath('~/MATdrr'));
flow=0;fhigh=125;dt=0.004;N=3;Niter=10;mode=1;verb=1;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d3=drr3drecon(d0,mask,flow,fhigh,dt,N,5,Niter,eps,verb,mode,a);
figure;imagesc([dc,dn,d0,d3]);%colormap(seis);
dl_snr(dc,d3)

%% USING DRR as initial model
param=struct('T',3,'niter',10,'mode',1,'K',64);
param.d0=d3;
mode=1;l1=6;l2=4;l3=1;s1=2;s2=2;s3=1;perc=2;Niter=5; %
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing for noisy data
% a=ones(Niter,1);
d4=dl_sgk_recon(d0,mask,mode,[l1,l2,l3],[s1,s2,s3],perc,Niter,a,param);
figure;imagesc([dc,dn,d0,d4]);%colormap(seis);
dl_snr(dc,d4) %






























