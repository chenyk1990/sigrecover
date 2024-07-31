%% THIS IS A SCRIPT TO DEMONSTRATE RESIDUAL DICTIONARY LEARNING FOR SIGNAL LEAKAGE RETRIEVAL
clc;clear;close all;

addpath(genpath('./subroutines'));


%% from urban DAS
%% patch size l1*l2
l1=32;l2=32;s1=16;s2=16;
c1=32;c2=16;%size of the 1D cosine dictionary (if c2>c1, overcomplete)
%% DCT dictionary (dctmtx will generates orthogonal transform)
dct=zeros(c1,c2);
for k=0:1:c2-1
    V=cos([0:1:c1-1]'*k*pi/c2);
    if k>0
        V=V-mean(V);
    end
    dct(:,k+1)=V/norm(V);
end
DCT=kron(dct,dct);%2D DCT dictionary (64,256)
param=struct('T',3,'niter',10,'mode',1,'K',64,'D',DCT);
perc=7;

load data/wang_in.mat;
load data/wang_out.mat;
d=d_in;
d3=d_out;
figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
imagesc([d,d3,d-d3]);colormap(sr_seis);caxis([-0.01,0.01]);

%test clip
d1=sr_clip(d,-0.04,0.04);
figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
imagesc([d,d1,d-d1]);colormap(sr_seis);caxis([-0.01,0.01]);

%test bp
d2=sr_bandpass(d,0.004,0,40);
figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
imagesc([d,d2,d-d2]);colormap(sr_seis);caxis([-0.01,0.01]);


%% residual learning
% X=sr_patch(d,mode,l1,l2,s1,s2);
% [D,G]=sr_sgk(X,param);
% G2=sr_ompN(D,X,T);
% perc=1;
% G2=sr_pthresh(G2,'ph',perc);
% X2=D*G2;
% d3=sr_patch_inv(X2,mode,n1,n2,l1,l2,s1,s2);
% figure;imagesc([dc,d,d3,d-d3]);colormap(sr_seis);
% sr_snr(dc,d3)
[n1,n2]=size(d);
XX=sr_patch(d3,1,l1,l2,s1,s2);
% XXn=sr_patch(sr_clip(d-d3,-0.02,0.02),1,l1,l2,s1,s2);
XXn=sr_patch(sr_bandpass(sr_clip(d-d3,-0.04,0.04),0.004,0,40),1,l1,l2,s1,s2);
[DD,GG]=sr_sgk(XX,param);
Gn=sr_ompN(DD,XXn,3);
perc=1;
Gn=sr_pthresh(Gn,'ph',perc);
Xn=DD*Gn;
d33=sr_patch_inv(Xn,1,n1,n2,l1,l2,s1,s2);
d33=sr_mf(d33,5,1,1);
d33=sr_mf(d33,5,1,2);
d4=d3+d33;
% figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
% imagesc([d,d3,d-d3,d33,d-d3-d33]);colormap(sr_seis);caxis([-0.01,0.01]);
% figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
% imagesc([d,d3,d4,d-d3,d-d4]);colormap(sr_seis);caxis([-0.01,0.01]);


perc=2;
d33_1=sr_patch_inv(DD*sr_pthresh(sr_ompN(DD,XXn,3),'ph',perc),1,n1,n2,l1,l2,s1,s2);
d33_1=sr_mf(d33_1,5,1,1);
d33_1=sr_mf(d33_1,5,1,1);

perc=0.5;
d33_2=sr_patch_inv(DD*sr_pthresh(sr_ompN(DD,XXn,3),'ph',perc),1,n1,n2,l1,l2,s1,s2);
d33_2=sr_mf(d33_2,5,1,1);
d33_2=sr_mf(d33_2,5,1,1);

perc=0.1;
d33_3=sr_patch_inv(DD*sr_pthresh(sr_ompN(DD,XXn,3),'ph',perc),1,n1,n2,l1,l2,s1,s2);
d33_3=sr_mf(d33_3,5,1,1);
d33_3=sr_mf(d33_3,5,1,1);


dt=0.20;%s
dx=1;%m
t=[0:n1-1]*dt;
x=1:n2;
indt=100:151;indx=900:1000;
indt2=200:250;indx2=200:400;

figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
subplot(2,3,1);imagesc(x,t,d33_1);colormap(sr_seis);caxis([-0.01,0.01]);title('Perc=2');ylabel('Time (s)','Fontsize',10,'fontweight','bold');sr_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
subplot(2,3,2);imagesc(x,t,d33_2);colormap(sr_seis);caxis([-0.01,0.01]);title("Perc=0.5");sr_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
subplot(2,3,3);imagesc(x,t,d33_3);colormap(sr_seis);caxis([-0.01,0.01]);title('Perc=0.1');sr_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');

% indt=1750:2000;indx=450:550;
subplot(2,3,4);imagesc(x(indx),t(indt),d33_1(indt,indx));colormap(sr_seis);caxis([-0.01,0.01]);title('Perc=2 (zoomed)');xlabel('Channel','Fontsize',10,'fontweight','bold');ylabel('Time (s)','Fontsize',10,'fontweight','bold');sr_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);sr_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
subplot(2,3,5);imagesc(x(indx),t(indt),d33_2(indt,indx));colormap(sr_seis);caxis([-0.01,0.01]);title("Perc=0.5 (zoomed)");xlabel('Channel','Fontsize',10,'fontweight','bold');sr_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);sr_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
subplot(2,3,6);imagesc(x(indx),t(indt),d33_3(indt,indx));colormap(sr_seis);caxis([-0.01,0.01]);title('Perc=0.1 (zoomed)');xlabel('Channel','Fontsize',10,'fontweight','bold');sr_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);sr_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
print(gcf,'-dpng','-r300','wang_perc.png');
print(gcf,'-depsc','-r300','wang_perc.eps');


