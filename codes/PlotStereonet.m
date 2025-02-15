function stereoplot=PlotStereonet(Dec,Inc,kwargs)

arguments
    Dec;
    Inc;
    kwargs.prjtype = 'equalarea';
    kwargs.fishermean {mustBeBoolean(kwargs.fishermean)} = false;
    kwargs.polargrid {mustBeBoolean(kwargs.polargrid)} = false;
    kwargs.saveplot {mustBeBoolean(kwargs.saveplot)} = false;
    kwargs.datacolor = 'blue';
    kwargs.datamarker = 'o';
    kwargs.datasize = 75;
    kwargs.holdonplots {mustBeBoolean(kwargs.holdonplots)} = true;
    kwargs.displaymean {mustBeBoolean(kwargs.displaymean)} = false;
    kwargs.displayprjtype {mustBeBoolean(kwargs.displayprjtype)} = false;
end

project1=strcmpi(kwargs.prjtype,{'equalangle','angle'});
project2=strcmpi(kwargs.prjtype,{'equalarea','area'});

N=length(Dec);

anglediff=nan(N,1);
Dec_f=nan(N,1);
Inc_f=nan(N,1);

for i=1:N

    anglediff(i,1)=AngDiff(Dec(1),Inc(1),Dec(i),Inc(i));

    if anglediff(i,1)>90
        Dec_f(i,1)=mod(Dec(i)+180,360);
        Inc_f(i,1)=-Inc(i);
    else
        Dec_f(i,1)=Dec(i);
        Inc_f(i,1)=Inc(i);
    end

end

[DecM,IncM,A95M,~,~]=FisherMean(Dec_f,Inc_f);
[Inc_sc,Dec_sc]=scircle1(IncM,DecM,A95M);

id=(Inc>=0);
idm=(IncM>=0);

if kwargs.holdonplots
    stereoplot=figure(1);
else
    stereoplot=figure;
end

set(stereoplot,'Position',[300,300,800,800]);

if sum(project1)==1

    projection='stereo';
    projorigin=[-90 0 0];
    
    if kwargs.displayprjtype
    plottitle='Equal-Angle Lower Hemisphere';
    end

end

if sum(project2)==1

    projection='eqdazim';
    projorigin=[-90 180 0];

    if kwargs.displayprjtype
    plottitle='Equal-Area Lower Hemisphere';
    end

end

if kwargs.polargrid

    axesm(projection,'frame','on','flinewidth',2,'ffacecolor','white','grid','on','glinestyle','-','glinewidth',0.5,'mlinelocation',10,'plinelocation',10,'origin',projorigin,'maplatlimit',[-90 0]);
    scatterm(-Inc(id),Dec(id),kwargs.datasize,kwargs.datamarker,'MarkerEdgeColor','black','MarkerFaceColor',kwargs.datacolor,'LineWidth',1.5);
    scatterm(Inc(~id),Dec(~id),kwargs.datasize,kwargs.datamarker,'MarkerEdgeColor',kwargs.datacolor,'MarkerFaceColor','white','LineWidth',1.5);

    if kwargs.fishermean
        plotm(-abs(Inc_sc),Dec_sc,'--','LineWidth',2,'Color',kwargs.datacolor);
        scatterm(-IncM(idm),DecM(idm),kwargs.datasize*3,'s','MarkerEdgeColor','black','MarkerFaceColor',kwargs.datacolor,'LineWidth',1.5);
        scatterm(IncM(~idm),DecM(~idm),kwargs.datasize*3,'s','MarkerEdgeColor',kwargs.datacolor,'MarkerFaceColor','white','LineWidth',1.5);
    end
    
    if kwargs.displayprjtype
    title(plottitle,'FontSize',15);
    end

    if kwargs.displaymean
        subtitle(['Dec = ',num2str(DecM,'%.1f'),'°, Inc = ',num2str(IncM,'%.1f'),'°, a95 = ',num2str(A95M,'%.1f'),'°, N = ',num2str(N)],'FontSize',15);
    end

end

if ~kwargs.polargrid

    axesm(projection,'frame','off','grid','off','origin',projorigin,'maplatlimit',[-90 0]);

    if kwargs.fishermean
        plotm(-abs(Inc_sc),Dec_sc,'--','LineWidth',2,'Color',kwargs.datacolor);
        scatterm(-IncM(idm),DecM(idm),kwargs.datasize*3,'s','MarkerEdgeColor','black','MarkerFaceColor',kwargs.datacolor,'LineWidth',1.5);
        scatterm(IncM(~idm),DecM(~idm),kwargs.datasize*3,'s','MarkerEdgeColor',kwargs.datacolor,'MarkerFaceColor','white','LineWidth',1.5);
    end

    scatterm(-Inc(id),Dec(id),kwargs.datasize,kwargs.datamarker,'MarkerEdgeColor','black','MarkerFaceColor',kwargs.datacolor,'LineWidth',1.5);
    scatterm(Inc(~id),Dec(~id),kwargs.datasize,kwargs.datamarker,'MarkerEdgeColor',kwargs.datacolor,'MarkerFaceColor','white','LineWidth',1.5);

    axesm(projection,'frame','on','flinewidth',2,'grid','on','glinestyle','-','glinewidth',0.5,'mlinelocation',10,'plinelocation',10,'origin',[0 0 0],'maplatlimit',[-90 90]);

    plotorder=get(gca,'Children');
    set(gca,'Children',flipud(plotorder));
    
    if kwargs.displayprjtype
    title(plottitle,'FontSize',15);
    end
    
    if kwargs.displaymean
        subtitle(['Dec = ',num2str(DecM,'%.1f'),'°, Inc = ',num2str(IncM,'%.1f'),'°, a95 = ',num2str(A95M,'%.1f'),'°, N = ',num2str(N)],'FontSize',15);
    end

end

tightmap;

if kwargs.saveplot
    exportgraphics(stereoplot,[plottitle,'.pdf'],'ContentType','vector');
end

end