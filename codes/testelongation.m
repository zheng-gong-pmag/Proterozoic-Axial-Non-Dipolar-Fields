%%

cd '/Users/zhenggong/Desktop/myfunctions';
close all; clear; clc;


Slat=48; Slon=277;
Plat=47; Plon=58; K=100; N=2000;

% Slat=-90+180*rand(1,1); 
% Slon=360*rand(1,1);
% Plat=-90+180*rand(1,1); 
% Plon=360*rand(1,1);
% K=100; 
% N=2000;

[Plon_r,Plat_r]=RandFisherDirs(Plon,Plat,K,N);

[Dec_r,Inc_r,~]=Pole2Dir(Slat,Slon,Plat_r,Plon_r,1);

[Dec,Inc,~,~,~]=FisherMean(Dec_r,Inc_r);

Latitude_r=atand(0.5*tand(Inc_r));

%% G3 terms

fig1=figure(1);

t1=tiledlayout(2,7,'TileSpacing','tight');

G3_list=-0.3:0.1:0.3;

for i=1:length(G3_list)

    G3=G3_list(i);
    [Inc_c]=getInclinationG2G3(Latitude_r,0,G3);

    % statistics on directions
    [Dec_f,Inc_f,~,~,~]=FisherMean(Dec_r,Inc_c);
    dpars = DoKent(Dec_r,Inc_c,length(Dec_r),'f');
    [Dec_b,Inc_b,Eig]=getEigen(Dec_r,Inc_c);
    [~,rnk]=max(Eig);

    Dec_fisher(i,1)=Dec_f;
    Inc_fisher(i,1)=Inc_f;

    Dec_kent(i,1)=dpars.dec;
    Inc_kent(i,1)=dpars.inc;

    Dec_bing(i,1)=Dec_b(rnk);
    Inc_bing(i,1)=Inc_b(rnk);

    % statistics on poles
    [Plat_c,Plon_c,~,~,~,~]=Dir2Pole(Slat,Slon,Dec_r,Inc_c,1);

    [Plon_f,Plat_f,~,K_f,~]=FisherMean(Plon_c,Plat_c);
    ppars = DoKent(Plon_c,Plat_c,length(Plon_c),'f');
    [Plon_b,Plat_b,Eig]=getEigen(Plon_c,Plat_c);
    [~,rnk]=max(Eig);
    Eig=sort(Eig,'descend');

    Plon_fisher(i,1)=Plon_f;
    Plat_fisher(i,1)=Plat_f;
    K_fisher(i,1)=K_f;

    Plon_kent(i,1)=ppars.dec;
    Plat_kent(i,1)=ppars.inc;
    K_kent(i,1)=ppars.Kappa;
    E_kent(i,1)=ppars.Ecc;
    
    Plon_bing(i,1)=Plon_b(rnk);
    Plat_bing(i,1)=Plat_b(rnk);
    K_bing(i,1)=Eig(1);
    E_bing(i,1)=Eig(3)./Eig(2);

%     nexttile(i);
%     PlotStereonet(Dec_r,Inc_c,'size',30,'polargrid',1,'color','red');
%     title(['G3 = ',num2str(G3*100),'%;'],'FontSize',15);
%     nexttile(i+7);
%     PlotStereonet(Plon_c,Plat_c,'size',30,'polargrid',1);
%     title(['G3 = ',num2str(G3*100),'%;'],'FontSize',15);

end

set(fig1,'Position',[300,400,1000,400]);

%%
fig2=figure(2);
set(fig2,'Position',[300,75,1000,700]);

subplot(2,3,1)
plot(G3_list*100,Inc*ones(1,length(G3_list)),'-k','LineWidth',2);
hold on
plot(G3_list*100,Inc_fisher,'-or','LineWidth',2);
plot(G3_list*100,Inc_kent,'-ob','LineWidth',2);
plot(G3_list*100,Inc_bing,'-og','LineWidth',2);
xlabel('G3 (%; octupole term)');ylabel('Inclination (°)');
xlim([-30,30]);xticks(-30:10:30);
% ylim([0,45]);
legend('initial value','Fisher statistics','Kent statistics','Location','northwest');

subplot(2,3,2)
plot(G3_list*100,Plat*ones(1,length(G3_list)),'-k','LineWidth',2);
hold on
plot(G3_list*100,Plat_fisher,'-or','LineWidth',2);
plot(G3_list*100,Plat_kent,'-ob','LineWidth',2);
plot(G3_list*100,Plat_bing,'-og','LineWidth',2);
xlabel('G3 (%; octupole term)');ylabel('Pole Latitude (°)');
xlim([-30,30]);xticks(-30:10:30);
% ylim([30,60]);

subplot(2,3,3)
plot(G3_list*100,Plon*ones(1,length(G3_list)),'-k','LineWidth',2);
hold on
plot(G3_list*100,Plon_fisher,'-or','LineWidth',2);
plot(G3_list*100,Plon_kent,'-ob','LineWidth',2);
plot(G3_list*100,Plon_bing,'-og','LineWidth',2);
xlabel('G3 (%; octupole term)');ylabel('Pole Longitude (°)');
xlim([-30,30]);xticks(-30:10:30);
% ylim([45,75]);

subplot(2,3,4)
plot(G3_list*100,K*ones(1,length(G3_list)),'-k','LineWidth',2);
hold on
plot(G3_list*100,K_fisher,'-or','LineWidth',2);
plot(G3_list*100,K_kent,'-ob','LineWidth',2);
plot(G3_list*100,K_bing,'-og','LineWidth',2);
xlabel('G3 (%; octupole term)');ylabel('Precision parameter');
xlim([-30,30]);xticks(-30:10:30);

subplot(2,3,5)
plot(G3_list*100,E_kent,'-ob','LineWidth',2);
hold on
plot(G3_list*100,E_bing,'-og','LineWidth',2);
xlabel('G3 (%; octupole term)');ylabel('Eccentricity');
xlim([-30,30]);xticks(-30:10:30);
ylim([0,1]);

subplot(2,3,6)
plot(K_fisher,E_kent,'-or','LineWidth',2);
hold on
plot(K_kent,E_kent,'-ob','LineWidth',2);
plot(K_bing,E_bing,'-og','LineWidth',2);
xlabel('Precision parameter');ylabel('Eccentricity');
% xlim([0,120]);xticks(0:20:120);

set(findall(gcf,'-property','FontSize'),'FontSize',12,'FontWeight','bold');


%% G2 terms

figure;

t2=tiledlayout(2,7,'TileSpacing','tight');

G2_list=-0.15:0.05:0.15;

for i=1:length(G2_list)

    G2=G2_list(i);
    [Inc_c]=getInclinationG2G3(Latitude_r,G2,0);

    % statistics on directions
    [Dec_f,Inc_f,~,~,~]=FisherMean(Dec_r,Inc_c);
    dpars = DoKent(Dec_r,Inc_c,length(Dec_r),'f');

    Dec_fisher(i,1)=Dec_f;
    Inc_fisher(i,1)=Inc_f;

    Dec_kent(i,1)=dpars.dec;
    Inc_kent(i,1)=dpars.inc;

    % statistics on poles
    [Plat_c,Plon_c,~,~,~,~]=Dir2Pole(Slat,Slon,Dec_r,Inc_c,1);

    [Plon_f,Plat_f,~,K_f,~]=FisherMean(Plon_c,Plat_c);
    ppars = DoKent(Plon_c,Plat_c,length(Plon_c),'f');

    Plon_fisher(i,1)=Plon_f;
    Plat_fisher(i,1)=Plat_f;
    K_fisher(i,1)=K_f;

    Plon_kent(i,1)=ppars.dec;
    Plat_kent(i,1)=ppars.inc;
    K_kent(i,1)=ppars.Kappa;
    E_kent(i,1)=ppars.Ecc;

%     nexttile(i);
%     PlotStereonet(Dec_r,Inc_c,'size',30,'polargrid',1,'color','red');
%     title(['G2 = ',num2str(G2*100),'%;'],'FontSize',15);
%     nexttile(i+7);
%     PlotStereonet(Plon_c,Plat_c,'size',30,'polargrid',1);
%     title(['G2 = ',num2str(G2*100),'%;'],'FontSize',15);

end

set(gcf,'Position',[300,400,1000,400]);

%%
fig4=figure(4);
set(fig4,'Position',[300,75,1000,700]);

subplot(2,3,1)
plot(G2_list*100,Inc*ones(1,length(G2_list)),'-k','LineWidth',2);
hold on
plot(G2_list*100,Inc_fisher,'-or','LineWidth',2);
plot(G2_list*100,Inc_kent,'-ob','LineWidth',2);
xlabel('G2 (%; quadrupole term)');ylabel('Inclination (°)');
xlim([-15,15]);xticks(-15:5:15);
% ylim([0,45]);
legend('initial value','Fisher statistics','Kent statistics','Location','northwest');

subplot(2,3,2)
plot(G2_list*100,Plat*ones(1,length(G2_list)),'-k','LineWidth',2);
hold on
plot(G2_list*100,Plat_fisher,'-or','LineWidth',2);
plot(G2_list*100,Plat_kent,'-ob','LineWidth',2);
xlabel('G2 (%; quadrupole term)');ylabel('Pole Latitude (°)');
xlim([-15,15]);xticks(-15:5:15);
% ylim([30,60]);

subplot(2,3,3)
plot(G2_list*100,Plon*ones(1,length(G2_list)),'-k','LineWidth',2);
hold on
plot(G2_list*100,Plon_fisher,'-or','LineWidth',2);
plot(G2_list*100,Plon_kent,'-ob','LineWidth',2);
xlabel('G2 (%; quadrupole term)');ylabel('Pole Longitude (°)');
xlim([-15,15]);xticks(-15:5:15);
% ylim([45,75]);

subplot(2,3,4)
plot(G2_list*100,K*ones(1,length(G2_list)),'-k','LineWidth',2);
hold on
plot(G2_list*100,K_fisher,'-or','LineWidth',2);
plot(G2_list*100,K_kent,'-ob','LineWidth',2);
xlabel('G2 (%; quadrupole term)');ylabel('Precision parameter');
xlim([-15,15]);xticks(-15:5:15);

subplot(2,3,5)
plot(G2_list*100,E_kent,'-ob','LineWidth',2);
xlabel('G2 (%; quadrupole term)');ylabel('Eccentricity');
xlim([-15,15]);xticks(-15:5:15);
ylim([0,1]);

subplot(2,3,6)
plot(K_fisher,E_kent,'-or','LineWidth',2);
hold on
plot(K_kent,E_kent,'-ob','LineWidth',2);
xlabel('Precision parameter');ylabel('Eccentricity');
% xlim([0,120]);xticks(0:20:120);

set(findall(gcf,'-property','FontSize'),'FontSize',12,'FontWeight','bold');
%%
% exportgraphics(fig1,['/Users/zhenggong/Desktop/1.pdf'],'ContentType','vector');
% exportgraphics(fig2,['/Users/zhenggong/Desktop/2.pdf'],'ContentType','vector');
% exportgraphics(fig3,['/Users/zhenggong/Desktop/3.pdf'],'ContentType','vector');
% exportgraphics(fig4,['/Users/zhenggong/Desktop/4.pdf'],'ContentType','vector');