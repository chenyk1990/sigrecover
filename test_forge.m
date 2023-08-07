
clc;clear;close all;



% load('~/dasdenoising/mat_raw/eq-3.mat');
% dn=d1;
% load('~/dasdenoising/mat_bpsomffk/eq-3.mat');
% d1=d1;

load('./eq-3.mat');
dn=d1;
load('./eq-3-d1.mat');
d1=d1;
figure;yc_imagesc([dn,d1,dn-d1]);
d2=yc_clip(dn-d1,-50,50);
figure;yc_imagesc([dn,d1,dn-d1,d2,dn-d1-d2]);





d=dn;
d3=d1;

figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
imagesc([d,d3,d-d3]);colormap(seis);caxis([-50,50]);

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
DCT=yc_initD([l1,l2,l3],[c1,c2,c3]);
param=struct('T',3,'niter',10,'mode',1,'K',64,'D',DCT);
perc=7;

[n1,n2]=size(d);
XX=yc_patch(d3,1,l1,l2,s1,s2);
% XXn=yc_patch(yc_clip(d-d3,-0.02,0.02),1,l1,l2,s1,s2);
XXn=yc_patch(yc_bandpass(yc_clip(d-d3,-80,80),0.004,0,60),1,l1,l2,s1,s2);
[DD,GG]=yc_sgk(XX,param);
Gn=yc_ompN(DD,XXn,3);
perc=1;
Gn=yc_pthresh(Gn,'ph',perc);
Xn=DD*Gn;
d33=yc_patch_inv(Xn,1,n1,n2,l1,l2,s1,s2);
% d33=yc_mf(d33,5,1,1);
% d33=yc_mf(d33,5,1,2);
d4=d3+d33;
% figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
% imagesc([d,d3,d-d3,d33,d-d3-d33]);colormap(seis);caxis([-50,50]);
% figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
% imagesc([d,d3,d4,d-d3,d-d4]);colormap(seis);caxis([-50,50]);

dt=0.005;%s
dx=1;%trace
t=[0:n1-1]*dt;
x=1:n2;
indt=300:600;indx=1:300;
indt2=1200:1400;indx2=740:960;
figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
subplot(2,3,1);imagesc(x,t,d);colormap(seis);caxis([-50,50]);title('Raw');ylabel('Time (s)','Fontsize',10,'fontweight','bold');yc_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);yc_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');text(-200,-1.2,'(a)','color','k','Fontsize',20,'fontweight','bold');
subplot(2,3,2);imagesc(x,t,d3);colormap(seis);caxis([-50,50]);title("IDF's");yc_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);yc_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');text(-200,-1.2,'(b)','color','k','Fontsize',20,'fontweight','bold');
subplot(2,3,3);imagesc(x,t,d4);colormap(seis);caxis([-50,50]);title('New');yc_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);yc_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');text(-200,-1.2,'(c)','color','k','Fontsize',20,'fontweight','bold');
subplot(2,3,4);imagesc(x,t,d33);colormap(seis);caxis([-50,50]);title('Retrieved signal');xlabel('Channel','Fontsize',10,'fontweight','bold');ylabel('Time (s)','Fontsize',10,'fontweight','bold');yc_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);yc_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');text(-200,-1.2,'(d)','color','k','Fontsize',20,'fontweight','bold');
subplot(2,3,5);imagesc(x,t,d-d3);colormap(seis);caxis([-50,50]);title("IDF's Noise");xlabel('Channel','Fontsize',10,'fontweight','bold');yc_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);yc_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');text(-200,-1.2,'(e)','color','k','Fontsize',20,'fontweight','bold');
subplot(2,3,6);imagesc(x,t,d-d4);colormap(seis);caxis([-50,50]);title('New noise');xlabel('Channel','Fontsize',10,'fontweight','bold');yc_framebox(x(indx(1)),x(indx(end)),t(indt(1)),t(indt(end)),'g',2);yc_framebox(x(indx2(1)),x(indx2(end)),t(indt2(1)),t(indt2(end)),'r',2);set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');text(-200,-1.2,'(f)','color','k','Fontsize',20,'fontweight','bold');
print(gcf,'-dpng','-r300','forge.png');
print(gcf,'-depsc','-r300','forge.eps');

% t=1:n1;
% x=1:n2;
% indt=1800:2000;indx=200:300;
figure('units','normalized','Position',[0.2 0.4 0.8, 1.0],'color','w');
subplot(4,3,1);imagesc(x(indx),t(indt),d(indt,indx));colormap(seis);caxis([-50,50]);title('Raw');ylabel('Time (s)','Fontsize',10,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');text(-70,1.2,'(a)','color','k','Fontsize',30,'fontweight','bold');
subplot(4,3,2);imagesc(x(indx),t(indt),d3(indt,indx));colormap(seis);caxis([-50,50]);title("IDF's");set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');
subplot(4,3,3);imagesc(x(indx),t(indt),d4(indt,indx));colormap(seis);caxis([-50,50]);title('New');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');
subplot(4,3,4);imagesc(x(indx),t(indt),d33(indt,indx));colormap(seis);caxis([-50,50]);title('Retrieved signal');ylabel('Time (s)','Fontsize',10,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');
subplot(4,3,5);imagesc(x(indx),t(indt),d(indt,indx)-d3(indt,indx));colormap(seis);caxis([-50,50]);title("IDF's Noise");set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');
subplot(4,3,6);imagesc(x(indx),t(indt),d(indt,indx)-d4(indt,indx));colormap(seis);caxis([-50,50]);title('New noise');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','g','ycolor','g');
% print(gcf,'-dpng','-r300','forge_z1.png');
%%

% t=1:n1;
% x=1:n2;
% indt2=1750:2000;indx2=450:550;
% figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
subplot(4,3,7);imagesc(x(indx2),t(indt2),d(indt2,indx2));colormap(seis);caxis([-50,50]);title('Raw');ylabel('Time (s)','Fontsize',10,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');text(690,5.8,'(b)','color','k','Fontsize',30,'fontweight','bold');
subplot(4,3,8);imagesc(x(indx2),t(indt2),d3(indt2,indx2));colormap(seis);caxis([-50,50]);title("IDF's");set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');
subplot(4,3,9);imagesc(x(indx2),t(indt2),d4(indt2,indx2));colormap(seis);caxis([-50,50]);title('New');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');
subplot(4,3,10);imagesc(x(indx2),t(indt2),d33(indt2,indx2));colormap(seis);caxis([-50,50]);title('Retrieved signal');xlabel('Channel','Fontsize',10,'fontweight','bold');ylabel('Time (s)','Fontsize',10,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');
subplot(4,3,11);imagesc(x(indx2),t(indt2),d(indt2,indx2)-d3(indt2,indx2));colormap(seis);caxis([-50,50]);xlabel('Channel','Fontsize',10,'fontweight','bold');title("IDF's noise");set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');
subplot(4,3,12);imagesc(x(indx2),t(indt2),d(indt2,indx2)-d4(indt2,indx2));colormap(seis);caxis([-50,50]);xlabel('Channel','Fontsize',10,'fontweight','bold');title('New noise');set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold','xcolor','r','ycolor','r');
% print(gcf,'-dpng','-r300','forge_z2.png');
print(gcf,'-dpng','-r300','forge_z.png');
print(gcf,'-depsc','-r300','forge_z.eps');
% % 
% 
% 

%% plot the atoms
figure('units','normalized','Position',[0.2 0.4 0.53, 1]);
for ia=1:64
    subplot(8,8,ia);imagesc(linspace(-1,1,l1),linspace(-1,1,l1),reshape(DCT(:,ia),l1,l2));colormap(seis);caxis([-0.2,0.2]);
    set(gca,'Linewidth',1.5,'Fontsize',16);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
end
print(gcf,'-depsc','-r300','forge_atoms1.eps');

figure('units','normalized','Position',[0.2 0.4 0.53, 1]);
rand('state',2024);
inds=randperm(64);
DD2=DD(:,inds);
inds=[11,14,16,18,19,21,23,26,27,28,33,35,36,38,43,49,53,59];
for ia=1:64
    subplot(8,8,ia);imagesc(linspace(-1,1,l1),linspace(-1,1,l1),reshape(DD2(:,ia),l1,l2));colormap(seis);caxis([-0.2,0.2]);
    if ismember(ia,inds)
        yc_draw_circle(0,0,0.5,'g.');
    end
    set(gca,'Linewidth',1.5,'Fontsize',16);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
end
print(gcf,'-depsc','-r300','forge_atoms2.eps');

save forge_atom.mat DCT DD DD2




% 
% 
% T=6;
% niter=12;
% K=256;
% n1=2000;
% n2=960;
% n3=1;
% ll1=32;
% ll2=16;
% ll3=1;
% ss1=8;
% ss2=8;
% ss3=1;
% perc=0.5;
% ifrobust=1;
% 
% %% simultaneous denoising and reconstruction
% %% SGK
% param=struct('T',T,'niter',niter,'mode',1,'K',K);
% mode=1;
% l1=ll1;
% l2=ll2;
% l3=ll3;
% s1=ss1;
% s2=ss2;
% s3=ss3;
% perc=perc;
% 
% %% create initial DCT
% c1=l1;c2=l2;c3=l3;
% Dinit=yc_initD([l1,l2,l3],[c1,c2,c3]);
% Dinit=Dinit(:,1:K);
% param.D=Dinit;
% 
% % %%test yc_initD
% % l1=4;l2=4;l3=1;
% % c1=l1;c2=l2;c3=l3;
% % Dinit=yc_initD([l1,l2,l3],[c1,c2,c3]);
% % % figure;imagesc(Dinit);
% 
% opt=1;%opt=1: update atoms or opt=0 not update
% % sigma=1.0;
% % epsi=1.345*sigma;
% epsi=50;
% % perc=100;
% % niter=10;niter=5;
% niter=5;
% if ifrobust==0
%     if n3==1
%         XX=yc_patch(d1,mode,l1,l2,s1,s2);
%         XXn=yc_patch(dn-d1,mode,l1,l2,s1,s2);
%     else
%         XX=yc_patch3d(d1,mode,l1,l2,l3,s1,s2,s3);
%         XXn=yc_patch3d(dn-d1,mode,l1,l2,l3,s1,s2,s3);
%     end
%     
%     
%     tic
%     [DD,GG]=yc_sgk(XX,param);
%     toc
%     
%     % for reference
%     %     tic
%     %     [DDksvd,GGksvd]=yc_ksvd(XX,param);
%     %     toc
%     
%     Gn=yc_ompN(DD,XXn,T);
%     Gn=yc_pthresh(Gn,'ph',perc);
%     Xn=DD*Gn;
%     if n3==1
%         d11=yc_patch_inv(Xn,mode,n1,n2,l1,l2,s1,s2);
%     else
%         d11=yc_patch3d_inv(Xn,mode,n1,n2,n3,l1,l2,l3,s1,s2,s3);
%     end
%     d2=d1+d11;
%     % figure;imagesc([dc,dn,d1,dn-d1,d11,dn-d1-d11]);colormap(seis);
%     % figure;imagesc([dc,dn,d4,d-d4]);colormap(seis);
%     
% else
%     if opt==0
%         if n3==1
%             XX=yc_patch(d1,mode,l1,l2,s1,s2);
%         else
%             XX=yc_patch3d(d1,mode,l1,l2,l3,s1,s2,s3);
%         end
%     end
%     tic
%     m=0;
%     d2=d1;
%     for it=1:niter
%         
%         n=dn-d1-m;
%         randn('state',232043124+it);
% %         n(find(n>epsi | n<-epsi ))=0.01*randn(size(find(n>epsi| n<-epsi)));
%         n=yc_clip(n,-epsi,epsi);
%         p=m+0.45*n;
% %         p=m+n;
%         
%         if n3==1
%             if opt==1
%                 XX=yc_patch(d2,mode,l1,l2,s1,s2);
%             end
%             XXn=yc_patch(p,mode,l1,l2,s1,s2);
%         else
%             if opt==1
%                 XX=yc_patch3d(d2,mode,l1,l2,l3,s1,s2,s3);
%             end
%             XXn=yc_patch3d(p,mode,l1,l2,l3,s1,s2,s3);
%         end
%         
%         
%         [DD,GG]=yc_sgk(XX,param);
%         
%         Gn=yc_ompN(DD,XXn,T);
%         Gn=yc_pthresh(Gn,'ph',perc);
%         Xn=DD*Gn;
%         if n3==1
%             d11=yc_patch_inv(Xn,mode,n1,n2,l1,l2,s1,s2);
%         else
%             d11=yc_patch3d_inv(Xn,mode,n1,n2,n3,l1,l2,l3,s1,s2,s3);
%         end
%         m=d11;
%         d2=d1+d11;
%         
%         fprintf("iteration = %d/%d \n",it,niter);
%         
%         figure(1);yc_imagesc([dn,d1,d11,d2,dn-d2]);pause(0.1);
%     end
%     toc
% end
% d2=reshape(d2,n1,n2*n3);
% 
% figure;yc_imagesc([dn,d1,d11,d2,dn-d2]);
% 
% 
% % rsf_create(denoised2,size(d2)');
% % rsf_write(d2,denoised2);
% 
% figure(1);yc_imagesc([dn,d1,n]);
% 
% 
% 
% 
% figure;yc_imagesc([dn,d1,dn-d1,d11,d2,dn-d2]);
% figure;yc_imagesc([dn,d1,dn-d1,d2,dn-d2]);
% 
% % 
% 
% % % % 
% nn=dn-d1;
% % dd=nn(200:600,1:200);
% dd=nn;
% % dd(find(nn>epsi | nn<-epsi ))=0.01*randn(size(find(nn>epsi| nn<-epsi)));
% % dd=yc_clip(dd,-epsi,epsi);
% [nn1,nn2]=size(dd);
% XXn=yc_patch(dd,mode,l1,l2,s1,s2);
% Gn=yc_ompN(DD,XXn,T);
% Gn=yc_pthresh(Gn,'ph',perc);
% Xn=DD*Gn;
% dd11=yc_patch_inv(Xn,mode,nn1,nn2,l1,l2,s1,s2);
% d3=d1+dd11;%non-robust version
% 
% 
% figure('units','normalized','Position',[0.2 0.4 0.55, 0.35],'color','w');
% yc_imagesc([dn,d1,dn-d1,d11,d2,dn-d2,dd11,d3,dn-d3]);
% print(gcf,'-depsc','-r300','figdas.eps');
% % save das.mat
% 
% % 
% % % figure;yc_imagesc([dd,dd11,dd-dd11]);
% % 
% % 
% % figure;yc_imagesc([nn,dd11,nn-dd11]);
% % 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
