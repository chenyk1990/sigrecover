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
figure('units','normalized','Position',[0.2 0.4 0.6, 1.0],'color','w');
subplot(6,3,1);imagesc(x,t,d);colormap(sr_seis);caxis([-0.01,0.01]);title('Raw');ylabel('Time (s)','Fontsize',10,'fontweight','bold');sr_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);sr_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');text(-220,-8.5,'(a)','color','k','Fontsize',20,'fontweight','bold');
subplot(6,3,2);imagesc(x,t,d3);colormap(sr_seis);caxis([-0.01,0.01]);title("Wang's");sr_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);sr_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');text(-220,-8.5,'(b)','color','k','Fontsize',20,'fontweight','bold');
subplot(6,3,3);imagesc(x,t,d4);colormap(sr_seis);caxis([-0.01,0.01]);title('New');sr_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);sr_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');text(-220,-8.5,'(c)','color','k','Fontsize',20,'fontweight','bold');
subplot(6,3,4);imagesc(x,t,d33);colormap(sr_seis);caxis([-0.01,0.01]);title('Retrieved signal');xlabel('Channel','Fontsize',10,'fontweight','bold');ylabel('Time (s)','Fontsize',10,'fontweight','bold');sr_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);sr_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');text(-220,-8.5,'(d)','color','k','Fontsize',20,'fontweight','bold');
subplot(6,3,5);imagesc(x,t,d-d3);colormap(sr_seis);caxis([-0.01,0.01]);title("Wang's Noise");xlabel('Channel','Fontsize',10,'fontweight','bold');sr_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);sr_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');text(-220,-8.5,'(e)','color','k','Fontsize',20,'fontweight','bold');
subplot(6,3,6);imagesc(x,t,d-d4);colormap(sr_seis);caxis([-0.01,0.01]);title('New noise');xlabel('Channel','Fontsize',10,'fontweight','bold');sr_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);sr_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');text(-220,-8.5,'(f)','color','k','Fontsize',20,'fontweight','bold');


% t=1:n1;
% x=1:n2;
% indt=1800:2000;indx=200:300;
% figure('units','normalized','Position',[0.2 0.4 0.8, 1.0],'color','w');
subplot(6,3,6+1);imagesc(x(indx),t(indt),d(indt,indx));colormap(sr_seis);caxis([-0.01,0.01]);title('Raw');ylabel('Time (s)','Fontsize',10,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');text(880,18,'(g)','color','k','Fontsize',20,'fontweight','bold');
subplot(6,3,6+2);imagesc(x(indx),t(indt),d3(indt,indx));colormap(sr_seis);caxis([-0.01,0.01]);title("Wang's");set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');
subplot(6,3,6+3);imagesc(x(indx),t(indt),d4(indt,indx));colormap(sr_seis);caxis([-0.01,0.01]);title('New');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');
subplot(6,3,6+4);imagesc(x(indx),t(indt),d33(indt,indx));colormap(sr_seis);caxis([-0.01,0.01]);title('Retrieved signal');ylabel('Time (s)','Fontsize',10,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');
subplot(6,3,6+5);imagesc(x(indx),t(indt),d(indt,indx)-d3(indt,indx));colormap(sr_seis);caxis([-0.01,0.01]);title("Wang's Noise");set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');
subplot(6,3,6+6);imagesc(x(indx),t(indt),d(indt,indx)-d4(indt,indx));colormap(sr_seis);caxis([-0.01,0.01]);title('New noise');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');
% print(gcf,'-dpng','-r300','wang_z1.png');
%%

% t=1:n1;
% x=1:n2;
% indt2=1750:2000;indx2=450:550;
% figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
subplot(6,3,6+7);imagesc(x(indx2),t(indt2),d(indt2,indx2));colormap(sr_seis);caxis([-0.01,0.01]);title('Raw');ylabel('Time (s)','Fontsize',10,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');text(160,38,'(h)','color','k','Fontsize',20,'fontweight','bold');
subplot(6,3,6+8);imagesc(x(indx2),t(indt2),d3(indt2,indx2));colormap(sr_seis);caxis([-0.01,0.01]);title("Wang's");set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');
subplot(6,3,6+9);imagesc(x(indx2),t(indt2),d4(indt2,indx2));colormap(sr_seis);caxis([-0.01,0.01]);title('New');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');
subplot(6,3,6+10);imagesc(x(indx2),t(indt2),d33(indt2,indx2));colormap(sr_seis);caxis([-0.01,0.01]);title('Retrieved signal');xlabel('Channel','Fontsize',10,'fontweight','bold');ylabel('Time (s)','Fontsize',10,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');
subplot(6,3,6+11);imagesc(x(indx2),t(indt2),d(indt2,indx2)-d3(indt2,indx2));colormap(sr_seis);caxis([-0.01,0.01]);xlabel('Channel','Fontsize',10,'fontweight','bold');title("Wang's noise");set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');
subplot(6,3,6+12);imagesc(x(indx2),t(indt2),d(indt2,indx2)-d4(indt2,indx2));colormap(sr_seis);caxis([-0.01,0.01]);xlabel('Channel','Fontsize',10,'fontweight','bold');title('New noise');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');
print(gcf,'-dpng','-r300','wang.png');
print(gcf,'-depsc','-r300','wang.eps');


figure('units','normalized','Position',[0.2 0.4 0.53, 1]);
for ia=1:64
    subplot(8,8,ia);imagesc(linspace(-1,1,l1),linspace(-1,1,l1),reshape(DCT(:,ia),l1,l2));colormap(sr_seis);caxis([-0.2,0.2]);
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
    subplot(8,8,ia);imagesc(linspace(-1,1,l1),linspace(-1,1,l1),reshape(DD2(:,ia),l1,l2));colormap(sr_seis);caxis([-0.2,0.2]);
    if ismember(ia,inds)
        sr_draw_circle(0,0,0.5,'g.');
    end
    set(gca,'Linewidth',1.5,'Fontsize',16);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
end
print(gcf,'-depsc','-r300','wang_atoms2.eps');
print(gcf,'-dpng','-r300','wang_atoms2.png');

save wang_atom.mat DCT DD DD2




%% plot patches
figure;imagesc(XXn(:,:));colormap(sr_seis);caxis([-0.01,0.01]);
figure;imagesc(Xn(:,:));colormap(sr_seis);caxis([-0.01,0.01]);



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
    subplot(8,8,ia-ntmp);imagesc(linspace(-1,1,l1),linspace(-1,1,l1),reshape(XXn(:,ia),l1,l1));colormap(sr_seis);caxis([-.01,.01]);
    if ismember(ia,inds)
        sr_draw_circle(0,0,0.5,'w');
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
    subplot(8,8,ia-ntmp);imagesc(linspace(-1,1,l1),linspace(-1,1,l1),reshape(Xn(:,ia),l1,l1));colormap(sr_seis);caxis([-.01,.01]);
    if ismember(ia,inds)
        sr_draw_circle(0,0,0.5,'w');
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
    subplot(8,8,ia-ntmp);imagesc(linspace(-1,1,l1),linspace(-1,1,l1),reshape(XXn(:,ia)-Xn(:,ia),l1,l1));colormap(sr_seis);caxis([-.01,.01]);
    if ismember(ia,inds)
        sr_draw_circle(0,0,0.5,'w');
    end
    set(gca,'Linewidth',1.5,'Fontsize',16);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
end
% print(gcf,'-depsc','-r300','wang_patch2.eps');


figure('units','normalized','Position',[0.2 0.4 0.53, 1]);
for ia=[1:16]
    subplot(4,4,ia);imagesc(linspace(-1,1,l1),linspace(-1,1,l1),reshape(XXn2(:,ia),l1,l1));colormap(sr_seis);caxis([-.02,.02]);
    set(gca,'Linewidth',1.5,'Fontsize',16);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
end
sr_anno2;
print(gcf,'-depsc','-r300','wang_patch0.eps');

figure('units','normalized','Position',[0.2 0.4 0.53, 1]);
for ia=[1:16]
    subplot(4,4,ia);imagesc(linspace(-1,1,l1),linspace(-1,1,l1),reshape(Xn2(:,ia),l1,l1));colormap(sr_seis);caxis([-.02,.02]);
    set(gca,'Linewidth',1.5,'Fontsize',16);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
end
sr_anno2;
print(gcf,'-depsc','-r300','wang_patch1.eps');

figure('units','normalized','Position',[0.2 0.4 0.53, 1]);
for ia=[1:16]
    subplot(4,4,ia);imagesc(linspace(-1,1,l1),linspace(-1,1,l1),reshape(XXn2(:,ia)-Xn2(:,ia),l1,l1));colormap(sr_seis);caxis([-.02,.02]);
    set(gca,'Linewidth',1.5,'Fontsize',16);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
end
sr_anno2;
print(gcf,'-depsc','-r300','wang_patch2.eps');


