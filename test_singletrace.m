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

%chenyk.data/students_postdocs/wanghang/urbandas_example/d_in.mat ./
%chenyk.data/students_postdocs/wanghang/urbandas_example/d_out.mat ./

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

dt=0.20;%s
dx=1;%m
t=[0:n1-1]*dt;
x=1:n2;
indt=100:151;indx=900:1000;indt2=200:250;indx2=200:400;

%% for single-trace comparison
idx=400;

figure('units','normalized','Position',[0.2 0.4 0.6, 1.0],'color','w');
subplot(4,3,1);imagesc(x,t,d);colormap(sr_seis);caxis([-0.01,0.01]);title('Raw','Fontsize',15,'fontweight','bold');ylabel('Time (s)','Fontsize',15,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');text(-220,-8.5,'(a)','color','k','Fontsize',20,'fontweight','bold');hold on;plot(x(idx)*ones(1,length(t)),t,'k--','linewidth',2);
subplot(4,3,2);imagesc(x,t,d3);colormap(sr_seis);caxis([-0.01,0.01]);title("Wang's",'Fontsize',15,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');text(-220,-8.5,'(b)','color','k','Fontsize',20,'fontweight','bold');hold on;plot(x(idx)*ones(1,length(t)),t,'k--','linewidth',2);
subplot(4,3,3);imagesc(x,t,d4);colormap(sr_seis);caxis([-0.01,0.01]);title('New','Fontsize',15,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');text(-220,-8.5,'(c)','color','k','Fontsize',20,'fontweight','bold');hold on;plot(x(idx)*ones(1,length(t)),t,'k--','linewidth',2);
subplot(4,3,4);imagesc(x,t,d33);colormap(sr_seis);caxis([-0.01,0.01]);title('Retrieved signal','Fontsize',15,'fontweight','bold');xlabel('Channel','Fontsize',15,'fontweight','bold');ylabel('Time (s)','Fontsize',15,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');text(-220,-8.5,'(d)','color','k','Fontsize',20,'fontweight','bold');hold on;plot(x(idx)*ones(1,length(t)),t,'k--','linewidth',2);
subplot(4,3,5);imagesc(x,t,d-d3);colormap(sr_seis);caxis([-0.01,0.01]);title("Wang's Noise",'Fontsize',15,'fontweight','bold');xlabel('Channel','Fontsize',15,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');text(-220,-8.5,'(e)','color','k','Fontsize',20,'fontweight','bold');
subplot(4,3,6);imagesc(x,t,d-d4);colormap(sr_seis);caxis([-0.01,0.01]);title('New noise','Fontsize',15,'fontweight','bold');xlabel('Channel','Fontsize',15,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');text(-220,-8.5,'(f)','color','k','Fontsize',20,'fontweight','bold');


subplot(2,1,2);plot(t,d(:,idx),'linewidth',2);
hold on;
plot(t,d3(:,idx),'linewidth',2);
plot(t,d4(:,idx),'linewidth',2);
plot(t,d33(:,idx),'linewidth',2);
legend('Raw',"Wang's",'New','Retreived','location','best','Fontsize',15,'Fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');
text(-220,0.07,'(g)','color','k','Fontsize',20,'fontweight','bold');
xlabel('Time (s)','Fontsize',15,'fontweight','bold');
ylabel('Amplitude','Fontsize',15,'fontweight','bold');
title('Single-trace comparison','Fontsize',15,'fontweight','bold');

print(gcf,'-dpng','-r300','wang_trace.png');


