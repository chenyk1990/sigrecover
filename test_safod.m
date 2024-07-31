
clc;clear;close all;


addpath(genpath('./subroutines'));

load('data/safod.mat');
dn=raw*10^10;dn=dn';
d1=d1*10^10;d1=d1';
figure;sr_imagesc([dn,d1,dn-d1]);
d2=sr_clip(dn-d1,-0.5,0.5);
figure;sr_imagesc([dn,d1,dn-d1,d2,dn-d1-d2]);

d=dn;
d3=d1;

figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
imagesc([d,d3,d-d3]);colormap(sr_seis);caxis([-0.5,0.5]);

%%
%% patch size l1*l2
l1=32;l2=16;l3=1;s1=16;s2=8;
l1=32;l2=16;l3=1;s1=8;s2=8;
c1=32;c2=16;c3=1;%size of the 1D cosine dictionary (if c2>c1, overcomplete)
% l1=48;l2=48;s1=24;s2=24;
% c1=48;c2=24;%size of the 1D cosine dictionary (if c2>c1, overcomplete)
%% DCT dictionary (dctmtx will generates orthogonal transform)
dct=zeros(c1,c2);
for k=0:1:c2-1
    V=cos([0:1:c1-1]'*k*pi/c2);
    if k>0
        V=V-mean(V);
    end
    dct(:,k+1)=V/norm(V);
end
% DCT=kron(dct,dct);%2D DCT dictionary (64,256)
DCT=sr_initD([l1,l2,l3],[c1,c2,c3]);
param=struct('T',3,'niter',10,'mode',1,'K',64,'D',DCT);
perc=7;

[n1,n2]=size(d);
XX=sr_patch(d3,1,l1,l2,s1,s2);
% XXn=sr_patch(sr_clip(d-d3,-0.02,0.02),1,l1,l2,s1,s2);
XXn=sr_patch(sr_bandpass(sr_clip(d-d3,-0.5,0.5),0.004,0,60),1,l1,l2,s1,s2);
[DD,GG]=sr_sgk(XX,param);
Gn=sr_ompN(DD,XXn,3);
perc=1;
Gn=sr_pthresh(Gn,'ph',perc);
Xn=DD*Gn;
d33=sr_patch_inv(Xn,1,n1,n2,l1,l2,s1,s2);
% d33=sr_mf(d33,5,1,1);
% d33=sr_mf(d33,5,1,2);
d4=d3+d33;
% figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
% imagesc([d,d3,d-d3,d33,d-d3-d33]);colormap(sr_seis);caxis([-0.5,0.5]);
% figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
% imagesc([d,d3,d4,d-d3,d-d4]);colormap(sr_seis);caxis([-0.5,0.5]);

dt=0.005;%s
dx=1;%trace
t=[0:n1-1]*dt;
x=1:n2;
indt=300:600;indx=1:300;
indt2=600:800;indx2=700:800;
figure('units','normalized','Position',[0.2 0.4 0.6, 1.0],'color','w');
subplot(6,3,1);imagesc(x,t,d);colormap(sr_seis);caxis([-0.5,0.5]);title('Raw');ylabel('Time (s)','Fontsize',10,'fontweight','bold');sr_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);sr_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');text(-200,-1.2,'(a)','color','k','Fontsize',20,'fontweight','bold');
subplot(6,3,2);imagesc(x,t,d3);colormap(sr_seis);caxis([-0.5,0.5]);title("IDF's");sr_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);sr_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');text(-200,-1.2,'(b)','color','k','Fontsize',20,'fontweight','bold');
subplot(6,3,3);imagesc(x,t,d4);colormap(sr_seis);caxis([-0.5,0.5]);title('New');sr_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);sr_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');text(-200,-1.2,'(c)','color','k','Fontsize',20,'fontweight','bold');
subplot(6,3,4);imagesc(x,t,d33);colormap(sr_seis);caxis([-0.5,0.5]);title('Retrieved signal');xlabel('Channel','Fontsize',10,'fontweight','bold');ylabel('Time (s)','Fontsize',10,'fontweight','bold');sr_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);sr_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');text(-200,-1.2,'(d)','color','k','Fontsize',20,'fontweight','bold');
subplot(6,3,5);imagesc(x,t,d-d3);colormap(sr_seis);caxis([-0.5,0.5]);title("IDF's Noise");xlabel('Channel','Fontsize',10,'fontweight','bold');sr_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);sr_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');text(-200,-1.2,'(e)','color','k','Fontsize',20,'fontweight','bold');
subplot(6,3,6);imagesc(x,t,d-d4);colormap(sr_seis);caxis([-0.5,0.5]);title('New noise');xlabel('Channel','Fontsize',10,'fontweight','bold');sr_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);sr_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');text(-200,-1.2,'(f)','color','k','Fontsize',20,'fontweight','bold');
% print(gcf,'-dpng','-r300','forge.png');


% t=1:n1;
% x=1:n2;
% indt=1800:2000;indx=200:300;
% figure('units','normalized','Position',[0.2 0.4 0.8, 1.0],'color','w');
subplot(6,3,6+1);imagesc(x(indx),t(indt),d(indt,indx));colormap(sr_seis);caxis([-0.5,0.5]);title('Raw');ylabel('Time (s)','Fontsize',10,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');text(-70,1.2,'(g)','color','k','Fontsize',20,'fontweight','bold');
subplot(6,3,6+2);imagesc(x(indx),t(indt),d3(indt,indx));colormap(sr_seis);caxis([-0.5,0.5]);title("IDF's");set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');
subplot(6,3,6+3);imagesc(x(indx),t(indt),d4(indt,indx));colormap(sr_seis);caxis([-0.5,0.5]);title('New');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');
subplot(6,3,6+4);imagesc(x(indx),t(indt),d33(indt,indx));colormap(sr_seis);caxis([-0.5,0.5]);title('Retrieved signal');ylabel('Time (s)','Fontsize',10,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');
subplot(6,3,6+5);imagesc(x(indx),t(indt),d(indt,indx)-d3(indt,indx));colormap(sr_seis);caxis([-0.5,0.5]);title("IDF's Noise");set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');
subplot(6,3,6+6);imagesc(x(indx),t(indt),d(indt,indx)-d4(indt,indx));colormap(sr_seis);caxis([-0.5,0.5]);title('New noise');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');


subplot(6,3,6+7);imagesc(x(indx2),t(indt2),d(indt2,indx2));colormap(sr_seis);caxis([-0.5,0.5]);title('Raw');ylabel('Time (s)','Fontsize',10,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');text(680,2.6,'(h)','color','k','Fontsize',20,'fontweight','bold');
subplot(6,3,6+8);imagesc(x(indx2),t(indt2),d3(indt2,indx2));colormap(sr_seis);caxis([-0.5,0.5]);title("IDF's");set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');
subplot(6,3,6+9);imagesc(x(indx2),t(indt2),d4(indt2,indx2));colormap(sr_seis);caxis([-0.5,0.5]);title('New');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');
subplot(6,3,6+10);imagesc(x(indx2),t(indt2),d33(indt2,indx2));colormap(sr_seis);caxis([-0.5,0.5]);title('Retrieved signal');xlabel('Channel','Fontsize',10,'fontweight','bold');ylabel('Time (s)','Fontsize',10,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');
subplot(6,3,6+11);imagesc(x(indx2),t(indt2),d(indt2,indx2)-d3(indt2,indx2));colormap(sr_seis);caxis([-0.5,0.5]);xlabel('Channel','Fontsize',10,'fontweight','bold');title("IDF's noise");set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');
subplot(6,3,6+12);imagesc(x(indx2),t(indt2),d(indt2,indx2)-d4(indt2,indx2));colormap(sr_seis);caxis([-0.5,0.5]);xlabel('Channel','Fontsize',10,'fontweight','bold');title('New noise');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');

print(gcf,'-dpng','-r300','safod.png');
print(gcf,'-depsc','-r300','safod.eps');