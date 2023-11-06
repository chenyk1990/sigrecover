%% THIS IS A SCRIPT TO DEMONSTRATE RESIDUAL DICTIONARY LEARNING FOR SIGNAL LEAKAGE RETRIEVAL
clc;clear;close all;

addpath(genpath('/Users/chenyk/chenyk/matlibcyk'));
addpath(genpath('./subroutines'));
% 
% %% DEMO script for DL (K-SVD) based denoising
% ne=20;%number of events;
% dt = 4./1000;
% tmax = dt*127;
% h = [-160:10:150];
% tau=linspace(0.1,1.8,ne);
% v0=linspace(1800,3000,ne);%exact
% randn('state',201819202122);amp = randn(ne,1);%[1., -1.,1];
% f0 = 20;
% snr = 2000;%default snr=2
% L = 6;
% seed=201517;
% dc=hevents(dt,f0,tmax,h,tau,v0,amp,snr,L,seed);dc=yc_scale(dc,2);
% 
% [n1,n2]=size(dc);dt=0.004;
% t=[0:n1-1]*dt;x=1:n2;
% %%
% randn('state',201617);
% d=dc+0.1*randn(size(dc));
% [n1,n2]=size(d);
% 
% %% patch size l1*l2
% l1=8;l2=8;s1=4;s2=4;
% c1=8;c2=16;%size of the 1D cosine dictionary (if c2>c1, overcomplete)
% %% DCT dictionary (dctmtx will generates orthogonal transform)
% dct=zeros(c1,c2);
% for k=0:1:c2-1
%     V=cos([0:1:c1-1]'*k*pi/c2);
%     if k>0
%         V=V-mean(V);
%     end
%     dct(:,k+1)=V/norm(V);
% end
% DCT=kron(dct,dct);%2D DCT dictionary (64,256)
% 
% % %test yc_initD
% % DCT2=yc_initD([l1,l2,1],[c2,c2,1]);
% % figure;yc_imagesc([DCT,DCT2,DCT-DCT2]);
% 
% 
% %% Denoising by KSVD
% %param naming following Chen, 2017, GJI; Zhou et al., 2020
% param.T=3;      %sparsity level
% param.D=DCT;    %initial D
% param.niter=10; %number of K-SVD iterations to perform; default: 10
% param.mode=1;   %1: sparsity; 0: error
% %param.exact:   Exact K-SVD update or approximate
% param.K=64;     %number of atoms, dictionary size
% %for X=DG
% %size of X: MxN
% %size of D: MxK
% %size of G: KxN
% 
% %% Option 1: denoise only using the integrated function
% param=struct('T',3,'niter',10,'mode',1,'K',64,'D',DCT);
% mode=1;l1=8;l2=8;s1=4;s2=4;perc=7;
% % d1=yc_ksvd_denoise(d,mode,[l1,l2,1],[s1,s2,1],perc,param);
% % figure;imagesc([dc,d,d1,d-d1]);colormap(seis);
% % yc_snr(dc,d1)
% % 
% % %% Option 2: without initialized param.D
% % param=rmfield(param,'D');%param=struct('T',3,'niter',10,'mode',1,'K',64);
% % mode=1;l1=8;l2=8;s1=4;s2=4;perc=7;
% % d1=yc_ksvd_denoise(d,mode,[l1,l2,1],[s1,s2,1],perc,param);
% % figure;imagesc([dc,d,d1,d-d1]);colormap(seis);
% % yc_snr(dc,d1)
% % 
% %% compare performance of two dictionaries
% T=3;perc=7;
% 
% Xc=yc_patch(dc,mode,l1,l2,s1,s2);
% X=yc_patch(d,mode,l1,l2,s1,s2);
% 
% [Dc,Gc]=yc_sgk(Xc,param);
% G1=yc_ompN(Dc,X,T);
% G1=yc_pthresh(G1,'ph',perc);
% X1=Dc*G1;
% d1=yc_patch_inv(X1,mode,n1,n2,l1,l2,s1,s2);
% yc_snr(dc,d1)
% 
% 
% [D,G]=yc_sgk(X,param);
% G2=yc_ompN(D,X,T);
% G2=yc_pthresh(G2,'ph',perc);
% X2=D*G2;
% d2=yc_patch_inv(X2,mode,n1,n2,l1,l2,s1,s2);
% yc_snr(dc,d2)
% 
% 
% % X2=D*G2;
% % d1=yc_patch_inv(X2,mode,n1,n2,l1,l2,s1,s2);
% % figure;imagesc([dc,d,d1,d-d1]);colormap(seis);
% % yc_snr(dc,d1)
% 


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
% [D,G]=yc_sgk(X,param);
% G2=yc_ompN(D,X,T);
% perc=1;
% G2=yc_pthresh(G2,'ph',perc);
% X2=D*G2;
% d3=yc_patch_inv(X2,mode,n1,n2,l1,l2,s1,s2);
% figure;imagesc([dc,d,d3,d-d3]);colormap(seis);
% yc_snr(dc,d3)
[n1,n2]=size(d);
XX=yc_patch(d3,1,l1,l2,s1,s2);
% XXn=yc_patch(yc_clip(d-d3,-0.02,0.02),1,l1,l2,s1,s2);
XXn=yc_patch(yc_bandpass(yc_clip(d-d3,-0.04,0.04),0.004,0,40),1,l1,l2,s1,s2);
[DD,GG]=yc_sgk(XX,param);
Gn=yc_ompN(DD,XXn,3);
perc=1;
Gn=yc_pthresh(Gn,'ph',perc);
Xn=DD*Gn;
d33=yc_patch_inv(Xn,1,n1,n2,l1,l2,s1,s2);
d33=yc_mf(d33,5,1,1);
d33=yc_mf(d33,5,1,2);
d4=d3+d33;
% figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
% imagesc([d,d3,d-d3,d33,d-d3-d33]);colormap(seis);caxis([-0.01,0.01]);
% figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
% imagesc([d,d3,d4,d-d3,d-d4]);colormap(seis);caxis([-0.01,0.01]);

dt=0.20;%s
dx=1;%m
t=[0:n1-1]*dt;
x=1:n2;
indt=100:151;indx=900:1000;indt2=200:250;indx2=200:400;
figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
subplot(2,3,1);imagesc(x,t,d);colormap(seis);caxis([-0.01,0.01]);title('Raw');ylabel('Time (s)','Fontsize',10,'fontweight','bold');yc_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);yc_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');text(-200,-6,'(a)','color','k','Fontsize',20,'fontweight','bold');
subplot(2,3,2);imagesc(x,t,d3);colormap(seis);caxis([-0.01,0.01]);title("Wang's");yc_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);yc_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');text(-200,-6,'(b)','color','k','Fontsize',20,'fontweight','bold');
subplot(2,3,3);imagesc(x,t,d4);colormap(seis);caxis([-0.01,0.01]);title('New');yc_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);yc_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');text(-200,-6,'(c)','color','k','Fontsize',20,'fontweight','bold');
subplot(2,3,4);imagesc(x,t,d33);colormap(seis);caxis([-0.01,0.01]);title('Retrieved signal');xlabel('Channel','Fontsize',10,'fontweight','bold');ylabel('Time (s)','Fontsize',10,'fontweight','bold');yc_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);yc_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');text(-200,-6,'(d)','color','k','Fontsize',20,'fontweight','bold');
subplot(2,3,5);imagesc(x,t,d-d3);colormap(seis);caxis([-0.01,0.01]);title("Wang's Noise");xlabel('Channel','Fontsize',10,'fontweight','bold');yc_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);yc_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');text(-200,-6,'(e)','color','k','Fontsize',20,'fontweight','bold');
subplot(2,3,6);imagesc(x,t,d-d4);colormap(seis);caxis([-0.01,0.01]);title('New noise');xlabel('Channel','Fontsize',10,'fontweight','bold');yc_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);yc_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');text(-200,-6,'(f)','color','k','Fontsize',20,'fontweight','bold');
print(gcf,'-dpng','-r300','wang.png');
print(gcf,'-depsc','-r300','wang.eps');

% t=1:n1;
% x=1:n2;
% indt=1800:2000;indx=200:300;
figure('units','normalized','Position',[0.2 0.4 0.8, 1.0],'color','w');
subplot(4,3,1);imagesc(x(indx),t(indt),d(indt,indx));colormap(seis);caxis([-0.01,0.01]);title('Raw');ylabel('Time (s)','Fontsize',10,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');text(880,18,'(a)','color','k','Fontsize',30,'fontweight','bold');
subplot(4,3,2);imagesc(x(indx),t(indt),d3(indt,indx));colormap(seis);caxis([-0.01,0.01]);title("Wang's");set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');
subplot(4,3,3);imagesc(x(indx),t(indt),d4(indt,indx));colormap(seis);caxis([-0.01,0.01]);title('New');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');
subplot(4,3,4);imagesc(x(indx),t(indt),d33(indt,indx));colormap(seis);caxis([-0.01,0.01]);title('Retrieved signal');ylabel('Time (s)','Fontsize',10,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');
subplot(4,3,5);imagesc(x(indx),t(indt),d(indt,indx)-d3(indt,indx));colormap(seis);caxis([-0.01,0.01]);title("Wang's Noise");set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');
subplot(4,3,6);imagesc(x(indx),t(indt),d(indt,indx)-d4(indt,indx));colormap(seis);caxis([-0.01,0.01]);title('New noise');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');
% print(gcf,'-dpng','-r300','wang_z1.png');
%%

% t=1:n1;
% x=1:n2;
% indt2=1750:2000;indx2=450:550;
% figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
subplot(4,3,7);imagesc(x(indx2),t(indt2),d(indt2,indx2));colormap(seis);caxis([-0.01,0.01]);title('Raw');ylabel('Time (s)','Fontsize',10,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');text(160,38,'(b)','color','k','Fontsize',30,'fontweight','bold');
subplot(4,3,8);imagesc(x(indx2),t(indt2),d3(indt2,indx2));colormap(seis);caxis([-0.01,0.01]);title("Wang's");set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');
subplot(4,3,9);imagesc(x(indx2),t(indt2),d4(indt2,indx2));colormap(seis);caxis([-0.01,0.01]);title('New');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');
subplot(4,3,10);imagesc(x(indx2),t(indt2),d33(indt2,indx2));colormap(seis);caxis([-0.01,0.01]);title('Retrieved signal');xlabel('Channel','Fontsize',10,'fontweight','bold');ylabel('Time (s)','Fontsize',10,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');
subplot(4,3,11);imagesc(x(indx2),t(indt2),d(indt2,indx2)-d3(indt2,indx2));colormap(seis);caxis([-0.01,0.01]);xlabel('Channel','Fontsize',10,'fontweight','bold');title("Wang's noise");set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');
subplot(4,3,12);imagesc(x(indx2),t(indt2),d(indt2,indx2)-d4(indt2,indx2));colormap(seis);caxis([-0.01,0.01]);xlabel('Channel','Fontsize',10,'fontweight','bold');title('New noise');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');
% print(gcf,'-dpng','-r300','wang_z2.png');
print(gcf,'-dpng','-r300','wang_z.png');
print(gcf,'-depsc','-r300','wang_z.eps');


figure('units','normalized','Position',[0.2 0.4 0.53, 1]);
for ia=1:64
    subplot(8,8,ia);imagesc(linspace(-1,1,l1),linspace(-1,1,l1),reshape(DCT(:,ia),l1,l2));colormap(seis);caxis([-0.2,0.2]);
    set(gca,'Linewidth',1.5,'Fontsize',16);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
end
print(gcf,'-depsc','-r300','wang_atoms1.eps');

figure('units','normalized','Position',[0.2 0.4 0.53, 1]);
rand('state',2024);
inds=randperm(64);
DD2=DD(:,inds);
inds=[5,7,8,9,11,14,16,19,21,26,27,28,31,33,35,36,38,41,44,53,59,64];
for ia=1:64
    subplot(8,8,ia);imagesc(linspace(-1,1,l1),linspace(-1,1,l1),reshape(DD2(:,ia),l1,l2));colormap(seis);caxis([-0.2,0.2]);
    if ismember(ia,inds)
        yc_draw_circle(0,0,0.5,'g.');
    end
    set(gca,'Linewidth',1.5,'Fontsize',16);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
end
print(gcf,'-depsc','-r300','wang_atoms2.eps');

save wang_atom.mat DCT DD DD2




%% plot patches
figure;imagesc(XXn(:,:));colormap(seis);caxis([-0.01,0.01]);
figure;imagesc(Xn(:,:));colormap(seis);caxis([-0.01,0.01]);



XXn2=[];
Xn2=[];
figure('units','normalized','Position',[0.2 0.4 0.53, 1]);
rand('state',2024);
inds=randperm(64);
DD2=DD(:,inds);
inds=[2,4,14,16,19,21,24,27,29,30,33,35,38,42,44,45,47,52,53,59,60,61,63];
ntmp=400

inds2=[4,5,6,7,8,11,12,13,36,37,38,39,40,41,42,44];
for ia=[1:64]+ntmp
    subplot(8,8,ia-ntmp);imagesc(linspace(-1,1,l1),linspace(-1,1,l1),reshape(XXn(:,ia),l1,l1));colormap(seis);caxis([-.01,.01]);
    if ismember(ia,inds)
        yc_draw_circle(0,0,0.5,'w');
    end
    set(gca,'Linewidth',1.5,'Fontsize',16);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);

    if ismember(ia-ntmp,inds2);
        XXn2=[XXn2,XXn(:,ia)];
    end
end



% print(gcf,'-depsc','-r300','wang_patch0.eps');

figure('units','normalized','Position',[0.2 0.4 0.53, 1]);
rand('state',2024);
inds=randperm(64);
DD2=DD(:,inds);
inds=[2,4,14,16,19,21,24,27,29,30,33,35,38,42,44,45,47,52,53,59,60,61,63];
% ntmp=4000
for ia=[1:64]+ntmp
    subplot(8,8,ia-ntmp);imagesc(linspace(-1,1,l1),linspace(-1,1,l1),reshape(Xn(:,ia),l1,l1));colormap(seis);caxis([-.01,.01]);
    if ismember(ia,inds)
        yc_draw_circle(0,0,0.5,'w');
    end
    set(gca,'Linewidth',1.5,'Fontsize',16);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
    if ismember(ia-ntmp,inds2);
        Xn2=[Xn2,Xn(:,ia)];
    end
end
% print(gcf,'-depsc','-r300','wang_patch1.eps');

figure('units','normalized','Position',[0.2 0.4 0.53, 1]);
rand('state',2024);
inds=randperm(64);
DD2=DD(:,inds);
inds=[2,4,14,16,19,21,24,27,29,30,33,35,38,42,44,45,47,52,53,59,60,61,63];
% ntmp=4000
for ia=[1:64]+ntmp
    subplot(8,8,ia-ntmp);imagesc(linspace(-1,1,l1),linspace(-1,1,l1),reshape(XXn(:,ia)-Xn(:,ia),l1,l1));colormap(seis);caxis([-.01,.01]);
    if ismember(ia,inds)
        yc_draw_circle(0,0,0.5,'w');
    end
    set(gca,'Linewidth',1.5,'Fontsize',16);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
end
% print(gcf,'-depsc','-r300','wang_patch2.eps');


figure('units','normalized','Position',[0.2 0.4 0.53, 1]);
for ia=[1:16]
    subplot(4,4,ia);imagesc(linspace(-1,1,l1),linspace(-1,1,l1),reshape(XXn2(:,ia),l1,l1));colormap(seis);caxis([-.02,.02]);
    set(gca,'Linewidth',1.5,'Fontsize',16);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
end
anno2;
print(gcf,'-depsc','-r300','wang_patch0.eps');

figure('units','normalized','Position',[0.2 0.4 0.53, 1]);
for ia=[1:16]
    subplot(4,4,ia);imagesc(linspace(-1,1,l1),linspace(-1,1,l1),reshape(Xn2(:,ia),l1,l1));colormap(seis);caxis([-.02,.02]);
    set(gca,'Linewidth',1.5,'Fontsize',16);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
end
anno2;
print(gcf,'-depsc','-r300','wang_patch1.eps');

figure('units','normalized','Position',[0.2 0.4 0.53, 1]);
for ia=[1:16]
    subplot(4,4,ia);imagesc(linspace(-1,1,l1),linspace(-1,1,l1),reshape(XXn2(:,ia)-Xn2(:,ia),l1,l1));colormap(seis);caxis([-.02,.02]);
    set(gca,'Linewidth',1.5,'Fontsize',16);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
end
anno2;
print(gcf,'-depsc','-r300','wang_patch2.eps');


