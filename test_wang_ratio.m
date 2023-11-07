%% THIS IS A SCRIPT TO DEMONSTRATE RESIDUAL DICTIONARY LEARNING FOR SIGNAL LEAKAGE RETRIEVAL
clc;clear;close all;

% addpath(genpath('/Users/chenyk/chenyk/matlibcyk'));

% download https://github.com/chenyk1990/MATseisdl

addpath(genpath('./MATseisdl'));
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
imagesc([d,d3,d-d3]);colormap(seis);caxis([-0.01,0.01]);

%test clip
d1=yc_clip(d,-0.04,0.04);
figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
imagesc([d,d1,d-d1]);colormap(seis);caxis([-0.01,0.01]);

%test bp
d2=yc_bandpass(d,0.004,0,40);
figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
imagesc([d,d2,d-d2]);colormap(seis);caxis([-0.01,0.01]);


%% residual learning
% X=yc_patch(d,mode,l1,l2,s1,s2);
% [D,G]=dl_sgk(X,param);
% G2=yc_ompN(D,X,T);
% perc=1;
% G2=dl_pthresh(G2,'ph',perc);
% X2=D*G2;
% d3=yc_patch_inv(X2,mode,n1,n2,l1,l2,s1,s2);
% figure;imagesc([dc,d,d3,d-d3]);colormap(seis);
% yc_snr(dc,d3)
[n1,n2]=size(d);
XX=yc_patch(d3,1,l1,l2,s1,s2);
% XXn=yc_patch(yc_clip(d-d3,-0.02,0.02),1,l1,l2,s1,s2);
XXn=yc_patch(yc_bandpass(yc_clip(d-d3,-0.04,0.04),0.004,0,40),1,l1,l2,s1,s2);
[DD,GG]=dl_sgk(XX,param);
Gn=yc_ompN(DD,XXn,3);
perc=1;
Gn=dl_pthresh(Gn,'ph',perc);
Xn=DD*Gn;
d33=yc_patch_inv(Xn,1,n1,n2,l1,l2,s1,s2);
d33=yc_mf(d33,5,1,1);
d33=yc_mf(d33,5,1,2);
d4=d3+d33;
% figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
% imagesc([d,d3,d-d3,d33,d-d3-d33]);colormap(seis);caxis([-0.01,0.01]);
% figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
% imagesc([d,d3,d4,d-d3,d-d4]);colormap(seis);caxis([-0.01,0.01]);


perc=2;
d33_1=yc_patch_inv(DD*dl_pthresh(yc_ompN(DD,XXn,3),'ph',perc),1,n1,n2,l1,l2,s1,s2);
d33_1=yc_mf(d33_1,5,1,1);
d33_1=yc_mf(d33_1,5,1,1);

perc=0.5;
d33_2=yc_patch_inv(DD*dl_pthresh(yc_ompN(DD,XXn,3),'ph',perc),1,n1,n2,l1,l2,s1,s2);
d33_2=yc_mf(d33_2,5,1,1);
d33_2=yc_mf(d33_2,5,1,1);

perc=0.1;
d33_3=yc_patch_inv(DD*dl_pthresh(yc_ompN(DD,XXn,3),'ph',perc),1,n1,n2,l1,l2,s1,s2);
d33_3=yc_mf(d33_3,5,1,1);
d33_3=yc_mf(d33_3,5,1,1);


dt=0.20;%s
dx=1;%m
t=[0:n1-1]*dt;
x=1:n2;
indt=100:151;indx=900:1000;
indt2=200:250;indx2=200:400;

figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
subplot(2,3,1);imagesc(x,t,d33_1);colormap(seis);caxis([-0.01,0.01]);title('Perc=2');ylabel('Time (s)','Fontsize',10,'fontweight','bold');yc_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
subplot(2,3,2);imagesc(x,t,d33_2);colormap(seis);caxis([-0.01,0.01]);title("Perc=0.5");yc_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
subplot(2,3,3);imagesc(x,t,d33_3);colormap(seis);caxis([-0.01,0.01]);title('Perc=0.1');yc_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');

% indt=1750:2000;indx=450:550;
subplot(2,3,4);imagesc(x(indx),t(indt),d33_1(indt,indx));colormap(seis);caxis([-0.01,0.01]);title('Perc=2 (zoomed)');xlabel('Channel','Fontsize',10,'fontweight','bold');ylabel('Time (s)','Fontsize',10,'fontweight','bold');yc_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);yc_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
subplot(2,3,5);imagesc(x(indx),t(indt),d33_2(indt,indx));colormap(seis);caxis([-0.01,0.01]);title("Perc=0.5 (zoomed)");xlabel('Channel','Fontsize',10,'fontweight','bold');yc_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);yc_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
subplot(2,3,6);imagesc(x(indx),t(indt),d33_3(indt,indx));colormap(seis);caxis([-0.01,0.01]);title('Perc=0.1 (zoomed)');xlabel('Channel','Fontsize',10,'fontweight','bold');yc_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);yc_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
print(gcf,'-dpng','-r300','wang_perc.png');
print(gcf,'-depsc','-r300','wang_perc.eps');

% % t=1:n1;
% % x=1:n2;
% % indt=1800:2000;indx=200:300;
% figure('units','normalized','Position',[0.2 0.4 0.8, 1.0],'color','w');
% subplot(4,3,1);imagesc(x(indx),t(indt),d(indt,indx));colormap(seis);caxis([-0.01,0.01]);title('Raw');ylabel('Time (s)','Fontsize',10,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');
% subplot(4,3,2);imagesc(x(indx),t(indt),d3(indt,indx));colormap(seis);caxis([-0.01,0.01]);title("Wang's");set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');
% subplot(4,3,3);imagesc(x(indx),t(indt),d4(indt,indx));colormap(seis);caxis([-0.01,0.01]);title('New');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');
% subplot(4,3,4);imagesc(x(indx),t(indt),d33(indt,indx));colormap(seis);caxis([-0.01,0.01]);title('Retrieved signal');ylabel('Time (s)','Fontsize',10,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');
% subplot(4,3,5);imagesc(x(indx),t(indt),d(indt,indx)-d3(indt,indx));colormap(seis);caxis([-0.01,0.01]);title("Wang's Noise");set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');
% subplot(4,3,6);imagesc(x(indx),t(indt),d(indt,indx)-d4(indt,indx));colormap(seis);caxis([-0.01,0.01]);title('New noise');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');
% % print(gcf,'-dpng','-r300','wang_z1.png');
% %%
% 
% % t=1:n1;
% % x=1:n2;
% % indt2=1750:2000;indx2=450:550;
% % figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
% subplot(4,3,7);imagesc(x(indx2),t(indt2),d(indt2,indx2));colormap(seis);caxis([-0.01,0.01]);title('Raw');ylabel('Time (s)','Fontsize',10,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');
% subplot(4,3,8);imagesc(x(indx2),t(indt2),d3(indt2,indx2));colormap(seis);caxis([-0.01,0.01]);title("Wang's");set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');
% subplot(4,3,9);imagesc(x(indx2),t(indt2),d4(indt2,indx2));colormap(seis);caxis([-0.01,0.01]);title('New');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');
% subplot(4,3,10);imagesc(x(indx2),t(indt2),d33(indt2,indx2));colormap(seis);caxis([-0.01,0.01]);title('Retrieved signal');xlabel('Channel','Fontsize',10,'fontweight','bold');ylabel('Time (s)','Fontsize',10,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');
% subplot(4,3,11);imagesc(x(indx2),t(indt2),d(indt2,indx2)-d3(indt2,indx2));colormap(seis);caxis([-0.01,0.01]);xlabel('Channel','Fontsize',10,'fontweight','bold');title("Wang's noise");set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');
% subplot(4,3,12);imagesc(x(indx2),t(indt2),d(indt2,indx2)-d4(indt2,indx2));colormap(seis);caxis([-0.01,0.01]);xlabel('Channel','Fontsize',10,'fontweight','bold');title('New noise');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');
% % print(gcf,'-dpng','-r300','wang_z2.png');
% print(gcf,'-dpng','-r300','wang_z.png');
% print(gcf,'-depsc','-r300','wang_z.eps');
% 
% 


