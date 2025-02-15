%% load and prepare data

clear; clc;  close all;
tic

datadir='/Users/zhenggong/Desktop/pmag_data/';

plotlist=[9];

for xx=1:length(plotlist)

    close all;

    xxx=plotlist(xx);

    filenamelist=["(01) 5 Ma N C",          "(02) 5 Ma R C",    "(03) 5 Ma ALL",        "(04) 10 Ma N C",       "(05) 10 Ma R C",   "(06) 10 Ma ALL",...
        "(07) Karoo Ferrar ALL",  "(08) CAMP",        "(09) Franklin ALL",    "(10) Umkondo ALL",     "(11) CSDG",        "(12) Mackenzie",...
        "(13) 1270 Ma ALL",       "(14) Aland ALL",   "(15) Matachewan ALL",  "(16) Bushveld ALL",    "(17) Dharwar ALL"];

    filename=char(filenamelist(xxx));
    filename=filename(6:end);

    D=csvread([datadir,filename,'.csv'],1,1); %#ok<CSVRD>

    SiteAge=D(:,1);     SiteLat=D(:,2);     SiteLon=D(:,3);     SiteDecl=D(:,4);    SiteIncl=D(:,5);
    a95=D(:,6);         Weight=D(:,7);      EulerLat=D(:,8);    EulerLon=D(:,9);    EulerAng=D(:,10);

    cd('/Users/zhenggong/Desktop/codes');

    % filter paleomagnetic data
    idx=(Weight>0);
    % [SitePlat,~,~,~,~,~]=Dir2VGP(SiteLat,SiteLon,SiteDecl,SiteIncl,a95);
    idy=(abs(SiteIncl)<=90);
    idz=SiteAge<=5000;

    % % cutoff
    % [SitePlat,SitePlon,~,~,~,~]=Dir2VGP(SiteLat,SiteLon,SiteDecl,SiteIncl,a95);
    % [~,~,~,~,~,idy]=CutOff(SitePlat,SitePlon,'method','V');

    id=idx.*idy.*idz;  id=(id>0);

    SiteLat=SiteLat(id);    SiteLon=SiteLon(id);    SiteDecl=SiteDecl(id);    SiteIncl=SiteIncl(id);      a95=a95(id);
    Weight=Weight(id);      EulerLat=EulerLat(id);  EulerLon=EulerLon(id);    EulerAng=EulerAng(id);

    % rotate site locations to paleo site locations using Euler parameters obtained from rotating mean paleopole to the north pole
    [Mean_SiteLon,Mean_SiteLat,~,~,~]=FisherMean(SiteLon,SiteLat);
    [SitePlat,SitePlon,~,~,~,~]=Dir2VGP(SiteLat,SiteLon,SiteDecl,SiteIncl,a95);
    [Mean_SitePlon,Mean_SitePlat,~,~,~]=FisherMean(SitePlon,SitePlat);
    [ELat,ELon,EAng]=VGP2GeoPole(Mean_SitePlat,Mean_SitePlon,'north');
    [SiteLat_paleo,~]=EulerRot(SiteLat,SiteLon,ELat,ELon,EAng);

    %% figure 1 - plot site locations

    fig1=figure(1);

    set(fig1,'position',[0,0,600,600]);
    %     axesm('MapProjection','ortho','Origin',[Mean_SiteLat Mean_SiteLon 0],'Frame','on','FLineWidth',1,...
    %         'Grid','on','mlinelocation',30,'plinelocation',30,'mlinelimit',[-90 90])

    axesm('MapProjection','robinson','Origin',[0 0 0],'Frame','on','Grid','on')
    load coastlines;
    plotm(coastlat,coastlon,'-','color',[128/255,128/255,128/255],'LineWidth',0.5);
    scatterm(SiteLat,SiteLon,50,'o','MarkerEdgeColor','black','MarkerFaceColor','yellow','LineWidth',1.5);

    title1='Site locations';
    title2=[filename,'; No. of sites = ',num2str(length(SiteLat))];
    title({title1,title2},'FontSize',15,'FontWeight','bold','fontname','Helvetica');


    %% figure 2 - plot Fisher statistics of VGPs

    G2_Min=-0.15;  G2_Max=0.15;  G2_Step=0.025;
    G3_Min=-0.30;  G3_Max=0.30;  G3_Step=0.025;
    G2=G2_Min:G2_Step:G2_Max; G2_length=length(G2);
    G3=G3_Min:G3_Step:G3_Max; G3=-G3; G3_length=length(G3);

    N=length(Weight);
    Paleolat=nan(N,1);

    bsnum=2000;
    K=nan(G3_length,G2_length,bsnum);
    KentK=nan(G3_length,G2_length,bsnum);
    KentB=nan(G3_length,G2_length,bsnum);
    KentE=nan(G3_length,G2_length,bsnum);
    tau1=nan(G3_length,G2_length,bsnum);
    tau32=nan(G3_length,G2_length,bsnum);


    for j=1:G3_length

        for k=1:G2_length

            for i=1:N
                [Paleolat(i,1)]=Paleolatitudefast(SiteIncl(i,1),SiteLat_paleo(i,1),G2(k),G3(j));
            end

            [PoleLat,PoleLon,~] = Pole(Paleolat,SiteLat,SiteLon,SiteDecl);
            [PoleLat_r,PoleLon_r] = EulerRot(PoleLat,PoleLon,EulerLat,EulerLon,EulerAng);

            for n=1:bsnum

                rnd=randsample((1:N),N,true);

                [~,~,~,KM,~] = FisherMean(PoleLon_r(rnd),PoleLat_r(rnd));
                [~,~,~,~,~,~,~,~,~,~,Kappa,Beta] = KentMean(PoleLon_r(rnd),PoleLat_r(rnd),N,'f');
                [~,~,Eig]=getEigen(PoleLon_r(rnd),PoleLat_r(rnd));
                Eig=sort(Eig,'descend');

                K(j,k,n)=KM;
                KentK(j,k,n)=Kappa;
                KentB(j,k,n)=Beta;
                KentE(j,k,n)=2*Beta./Kappa;
                tau1(j,k,n)=Eig(1);
                tau32(j,k,n)=Eig(3)./Eig(2);

            end

            K_m(j,k)=median(reshape(K(j,k,:),[],1));

            KentK_m(j,k)=median(reshape(KentK(j,k,:),[],1));

            KentE_m(j,k)=median(reshape(KentE(j,k,:),[],1));

            tau1_m(j,k)=median(reshape(tau1(j,k,:),[],1));

            tau32_m(j,k)=median(reshape(tau32(j,k,:),[],1));

        end

        disp(j);

    end


    for j=1:G3_length

        for k=1:G2_length

            K_delta(j,k,:)=reshape(K(j,k,:),[],1)-reshape(K((G3_length+1)/2,(G2_length+1)/2,:),[],1);
            K_p(j,k)=zeroinCI2P(reshape(K_delta(j,k,:),[],1));

            KentK_delta(j,k,:)=reshape(KentK(j,k,:),[],1)-reshape(KentK((G3_length+1)/2,(G2_length+1)/2,:),[],1);
            KentK_p(j,k)=zeroinCI2P(reshape(KentK(j,k,:),[],1)-reshape(KentK((G3_length+1)/2,(G2_length+1)/2,:),[],1));

            KentE_delta(j,k,:)=reshape(KentE(j,k,:),[],1)-reshape(KentE((G3_length+1)/2,(G2_length+1)/2,:),[],1);
            KentE_p(j,k)=zeroinCI2P(reshape(KentE(j,k,:),[],1)-reshape(KentE((G3_length+1)/2,(G2_length+1)/2,:),[],1));

            tau1_delta(j,k,:)=reshape(tau1(j,k,:),[],1)-reshape(tau1((G3_length+1)/2,(G2_length+1)/2,:),[],1);
            tau1_p(j,k)=zeroinCI2P(reshape(tau1(j,k,:),[],1)-reshape(tau1((G3_length+1)/2,(G2_length+1)/2,:),[],1));

            tau32_delta(j,k,:)=reshape(tau32(j,k,:),[],1)-reshape(tau32((G3_length+1)/2,(G2_length+1)/2,:),[],1);
            tau32_p(j,k)=zeroinCI2P(reshape(tau32(j,k,:),[],1)-reshape(tau32((G3_length+1)/2,(G2_length+1)/2,:),[],1));

        end
    end

    %%
    fig2=figure(2);
    set(fig2,'Position',[20,250,1400,400]);
    title1='VGP Statistics';
    title2=[filename,'; No. of sites = ',num2str(length(SiteLat))];
    sgtitle({title1,title2},'FontSize',15,'FontWeight','bold','fontname','Helvetica');

    s1=subplot(151);
    pcolor(G2*100,G3*100,K_m);
    hold on
    shading interp
    [c,h]=contour(G2*100,G3*100,K_p,[0.05,0.1,0.25],"ShowText",true);
    clabel(c,h,'FontSize',15,'Color','w','fontname','Helvetica');
    set(h,'LineWidth',2,'LineColor','w','LineStyle','-');
    plot(0,0,'d','MarkerSize',10,'MarkerFaceColor','w');
    hold off
    subtitle('Fisher Precision (K)','FontSize',15,'fontname','Helvetica');
    xlabel('G2 (%; quadrupole term)');
    ylabel('G3 (%; octupole term)');
    colorbar('EastOutside');
    set(gca,'FontSize',12,'TickDir','out','fontname','Helvetica');
    daspect([1 1 1]);
    colormap(s1,flip(viridis));
    clim([min(K_m,[],'all') max(K_m,[],'all')]);

    s2=subplot(152);
    pcolor(G2*100,G3*100,tau1_m);
    hold on
    shading interp
    [c,h]=contour(G2*100,G3*100,tau1_p,[0.05,0.1,0.25],"ShowText",true);
    clabel(c,h,'FontSize',15,'Color','w','fontname','Helvetica');
    set(h,'LineWidth',2,'LineColor','w','LineStyle','-');
    plot(0,0,'d','MarkerSize',10,'MarkerFaceColor','w');
    hold off
    subtitle('Bingham Eigen (\tau\it1)','FontSize',15,'fontname','Helvetica');
    xlabel('G2 (%; quadrupole term)');
    ylabel('G3 (%; octupole term)');
    colorbar('EastOutside');
    set(gca,'FontSize',12,'TickDir','out','fontname','Helvetica');
    daspect([1 1 1]);
    colormap(s2,flip(viridis));
    clim([min(tau1_m,[],'all') max(tau1_m,[],'all')]);

    s3=subplot(153);
    pcolor(G2*100,G3*100,KentK_m);
    hold on
    shading interp
    [c,h]=contour(G2*100,G3*100,KentK_p,[0.05,0.1,0.25],"ShowText",true);
    clabel(c,h,'FontSize',15,'Color','w','fontname','Helvetica');
    set(h,'LineWidth',2,'LineColor','w','LineStyle','-');
    plot(0,0,'d','MarkerSize',10,'MarkerFaceColor','w');
    hold off
    subtitle('Kent Concentration (\kappa)','FontSize',15,'fontname','Helvetica');
    xlabel('G2 (%; quadrupole term)');
    ylabel('G3 (%; octupole term)');
    colorbar('EastOutside');
    set(gca,'FontSize',12,'TickDir','out','fontname','Helvetica');
    daspect([1 1 1]);
    colormap(s3,flip(viridis));
    clim([min(KentK_m,[],'all') max(KentK_m,[],'all')]);

    s4=subplot(154);
    pcolor(G2*100,G3*100,tau32_m);
    hold on
    shading interp
    [c,h]=contour(G2*100,G3*100,tau32_p,[0.05,0.1,0.25],"ShowText",true);
    clabel(c,h,'FontSize',15,'Color','w','fontname','Helvetica');
    set(h,'LineWidth',2,'LineColor','w','LineStyle','-');
    plot(0,0,'d','MarkerSize',10,'MarkerFaceColor','w');
    hold off
    subtitle('Bingham Eigen (\tau\it3/\tau\it2)','FontSize',15,'fontname','Helvetica');
    xlabel('G2 (%; quadrupole term)');
    ylabel('G3 (%; octupole term)');
    colorbar('EastOutside');
    set(gca,'FontSize',12,'TickDir','out','fontname','Helvetica');
    daspect([1 1 1]);
    colormap(s4,flip(viridis));
    clim([min(tau32_m,[],'all') max(tau32_m,[],'all')]);

    s5=subplot(155);
    pcolor(G2*100,G3*100,KentE_m);
    hold on
    shading interp
    [c,h]=contour(G2*100,G3*100,KentE_p,[0.05,0.1,0.25],"ShowText",true);
    clabel(c,h,'FontSize',15,'Color','w','fontname','Helvetica');
    set(h,'LineWidth',2,'LineColor','w','LineStyle','-');
    plot(0,0,'d','MarkerSize',10,'MarkerFaceColor','w');
    hold off
    subtitle('Kent Eccentricity (\it2\beta/\kappa)','FontSize',15,'fontname','Helvetica');
    xlabel('G2 (%; quadrupole term)');
    ylabel('G3 (%; octupole term)');
    colorbar('EastOutside');
    set(gca,'FontSize',12,'TickDir','out','fontname','Helvetica');
    daspect([1 1 1]);
    colormap(s5,viridis);
    clim([min(KentE_m,[],'all') max(KentE_m,[],'all')]);

    %% figure 3 - plot distribution of VGPs

    fig3=figure(3);
    set(fig3,'Position',[600,0,800,800]);
    title1='VGP Distribution';
    title2=[filename,'; No. of sites = ',num2str(length(SiteLat))];
    sgtitle({title1,title2},'FontSize',15,'FontWeight','bold','fontname','Helvetica');

    G2_Min=-0.15; G2_Max=0.15; G2_Step=0.15;
    G3_Min=-0.30; G3_Max=0.30; G3_Step=0.30;
    G2_length=length(G2_Min:G2_Step:G2_Max);
    G3_length=length(G3_Min:G3_Step:G3_Max);

    N=length(Weight);
    Paleolat=zeros(N,1);

    for j=1:G3_length
        G_3=G3_Max-G3_Step*(j-1);

        for k=1:G2_length
            G_2=G2_Min+G2_Step*(k-1);

            for i=1:N
                [Paleolat(i)]=Paleolatitudefast(SiteIncl(i,1),SiteLat_paleo(i,1),G_2,G_3);
            end

            [PoleLat,PoleLon,~] = Pole(Paleolat,SiteLat,SiteLon,SiteDecl);
            [PoleLat_r,PoleLon_r] = EulerRot(PoleLat,PoleLon,EulerLat,EulerLon,EulerAng);
            [MeanVGP_Lon,MeanVGP_Lat,Mean_A95,~,~] = FisherMean(PoleLon_r,PoleLat_r);
            [MeanVGP_Lat_sc,MeanVGP_Lon_sc]=scircle1(MeanVGP_Lat,MeanVGP_Lon,Mean_A95);

            subplot(G3_length,G2_length,(j-1)*G2_length+k)
            axesm('MapProjection','stereo','Origin',[(Mean_SiteLat+MeanVGP_Lat)/2 ,(Mean_SiteLon+MeanVGP_Lon)/2, 0],'Frame','on','FLineWidth',1,...
                'Grid','on','mlinelocation',30,'plinelocation',30,'mlinelimit',[-90 90])
            %                     axesm('MapProjection','ortho','Origin',[MeanVGP_Lat,MeanVGP_Lon,0],'Frame','on','FLineWidth',1,...
            %                         'Grid','on','mlinelocation',30,'plinelocation',30,'mlinelimit',[-90 90])
            %             axesm('MapProjection','mollweid','Origin',[0 0 0],'Frame','on','Grid','off')
            load coastlines;
            plotm(coastlat,coastlon,'-','color',[128/255,128/255,128/255],'LineWidth',0.5);
            scatterm(PoleLat_r,PoleLon_r,50,'o','MarkerEdgeColor','black','MarkerFaceColor','yellow','LineWidth',1);
            scatterm(MeanVGP_Lat,MeanVGP_Lon,75,'d','MarkerEdgeColor','black','MarkerFaceColor','red','LineWidth',1.5);
            scatterm(Mean_SiteLat,Mean_SiteLon,75,'d','MarkerEdgeColor','black','MarkerFaceColor','green','LineWidth',1.5);
            plotm(MeanVGP_Lat_sc,MeanVGP_Lon_sc,'-r','LineWidth',1);
            subtitle(['G3=',num2str(G_3.*100),'%; G2=',num2str(G_2.*100),'%'],'FontSize',15,'fontname','Helvetica');

        end
    end



    %% figure 4 - plot directions

    % fig4=figure(4);
    % set(gcf,'Position',get(0,'Screensize'));
    % t=tiledlayout(3,3,'TileSpacing','tight');
    % title1='Distribution of VGPs';
    % title2=[filename,'; No. of sites = ',num2str(length(SiteLat))];
    % sgtitle({title1,title2},'FontSize',15,'FontWeight','bold','fontname','Helvetica');
    %
    % G2_Min=-0.15; G2_Max=0.15; G2_Step=0.15;
    % G3_Min=-0.30; G3_Max=0.30; G3_Step=0.30;
    % G2_length=length(G2_Min:G2_Step:G2_Max);
    % G3_length=length(G3_Min:G3_Step:G3_Max);
    %
    % N=length(Weight);
    % Paleolat=zeros(N,1);
    %
    % for j=1:G3_length
    %     G_3=G3_Max-G3_Step*(j-1);
    %
    %     for k=1:G2_length
    %         G_2=G2_Min+G2_Step*(k-1);
    %
    %         for i=1:N
    %             [Paleolat(i)]=Paleolatitudefast(SiteIncl(i,1),SiteLat_paleo(i,1),G_2,G_3);
    %         end
    %
    %         [~,~,PaleoSiteIncl] = Pole(Paleolat,SiteLat,SiteLon,SiteDecl);
    %
    %         nexttile;
    %         PlotStereonet(SiteDecl,PaleoSiteIncl,'datasize',30,'polargrid',1);
    %     end
    %
    % end


    %% figure 5 - plot distributions

    % figure
    % i=1;
    % for j=1:6:25
    %
    %     for k=1:6:13
    %
    %         a=reshape(KentE_delta(j,k,:),[],1);
    %         a0=reshape(KentE(5,3,:),[],1);
    %
    %         subplot(5,3,i)
    %         histogram(a)
    %
    %         i=i+1;
    %     end
    % end



    %% export graphics
    % exportgraphics(fig1,'Fig1 Site and VGP locations.pdf','ContentType','vector');
    % exportgraphics(fig2,'Fig2 Fisher statistics of VGPs.pdf','ContentType','vector');
    % exportgraphics(fig3,'Fig3 Distribution of VGPs.pdf','ContentType','vector');
    %
    % exportgraphics(fig1,['/Users/zhenggong/Desktop/',num2str(xxx),' ',filename,'-a.pdf'],'ContentType','vector');
    % exportgraphics(fig2,['/Users/zhenggong/Desktop/',num2str(xxx),' ',filename,'-b.pdf'],'ContentType','vector');
    % exportgraphics(fig3,['/Users/zhenggong/Desktop/',num2str(xxx),' ',filename,'-c.pdf'],'ContentType','vector');

end
toc