close all; clc; clear;
% 
G2=-0.15:0.075:0.15;
G3=-(-0.3:0.15:0.30);

% G2=zeros(1,7);
% G3=0.15:0.05:0.45;


% G3=0.0;
% G2=-1.10000;

rad=pi/180;
phi=(0:1:180);
Lat=90-phi;
phi=phi.*rad;

GADIncl=atan(2.*cos(phi)./abs(sin(phi)))./rad;

%% Plot Latitude vs Inclination

figure
set(gcf,'Position',get(0,'Screensize'));

for j=1:length(G3)

    for k=1:length(G2)

        G_2=G2(k);G_3=G3(j);

        I_num=2.*cos(phi)+1.5.*G_2.*(3*cos(phi).^2-1)+2.*G_3.*(5.*cos(phi).^3-3.*cos(phi));
        I_den=sin(phi)+G_2.*(3.*sin(phi).*cos(phi))+1.5.*G_3.*(5.*sin(phi).*cos(phi).^2-sin(phi));
        
        Incl1=atan(I_num./I_den)./rad;
        Incl2=atan(I_num./abs(I_den))./rad;
        

        subplot(length(G3),length(G2),(j-1)*length(G2)+k)
        
        plot(Lat,GADIncl,'-k','LineWidth',1.5);
        hold on
        plot(Lat,Incl1,'-r','LineWidth',2.5);
        plot(Lat,Incl2,'-b','LineWidth',1.5);
        

        if G_2==-0.15
            ylabel('Inclination (°)');
        end

        if G_3==-0.3
            xlabel('Latitude (°)');
        end

        if j*k==1
            legend('GAD','non-GAD Eq.1','non-GAD Eq.2','Location','best');
        end

        xlim([-90 90]);ylim([-90 90]);
        xticks(-90:30:90); yticks(-90:30:90);
        subtitle(['G2=',num2str(G_2*100),'%; G3=',num2str(G_3*100),'%'],"FontWeight","bold");

    end

end

%% Plot delta(Latitude) vs Inclination

figure
set(gcf,'Position',get(0,'Screensize'));

for j=1:length(G3)

    for k=1:length(G2)

        G_2=G2(k);G_3=G3(j);

        I_num=2.*cos(phi)+1.5.*G_2.*(3*cos(phi).^2-1)+2.*G_3.*(5.*cos(phi).^3-3.*cos(phi));
        I_den=sin(phi)+G_2.*(3.*sin(phi).*cos(phi))+1.5.*G_3.*(5.*sin(phi).*cos(phi).^2-sin(phi));

        Incl2=atan(I_num./abs(I_den))./rad;
        % Incl=atan(I_num./I_den)./rad;

        L1=fit(Lat',GADIncl','linear');
        L2=fit(Lat',Incl2','linear');

        subplot(length(G3),length(G2),(j-1)*length(G2)+k)

        deltaLat=feval(L1,Incl2)-feval(L1,GADIncl);
%         plot(feval(L1,GADIncl),feval(L1,Incl),'-b','LineWidth',1.5);
                plot(deltaLat, Incl2,'-b','LineWidth',1.5);
        hold on
        plot([-90,90],[-90,90],'-r','LineWidth',1.);

        if G_2==-0.15
            ylabel('Inclination (°)');
        end

        if G_3==-0.3
            xlabel('\DeltaLatitude (°)');
        end

        %         if j*k==1
        %             legend('GAD','non-GAD','Location','best');
        %         end

        xlim([-45 45]);ylim([-90 90]);
        xticks(-90:10:90); yticks(-90:30:90);
        subtitle(['G2=',num2str(G_2*100),'%; G3=',num2str(G_3*100),'%'],"FontWeight","bold");

    end

end

%%
% exportgraphics(gcf,'/Users/zhenggong/Desktop/Inc vs Lat.pdf','ContentType','vector');
