%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% The following algorithm allows forward modeling of beheaded stream profiles for a %%
%%%%%% range of area loss scenarios. The underlying theory is presented in a manuscript %%%
%%%%%% which is currently under review for publication in Journal of Geophysical %%%%%%%%%%
%%%%%% Research - Earth Surface: "A new method to restore tectonically %%%%%%%%%%%%%%%%%%%%
%%%%%% beheaded valleys" by Adrien Moulin, Matthieu Ribot, and Sigurjón Jónsson. Under %%%%
%%%%%% its present form, the algorithm performs the inversion of area loss and tectonic %%%
%%%%%% uplift recorded by 2 beheaded streams across the Wadi-al-Akhdar graben %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In the following T, ST3A, ST3B, and ST3C refer to the stream sections defined in Moulin et al. in review
% confluence ST3A-B and confluence ST3C refer to the two branching points defined in the same article


%% Extract relevant information from the TopoToolBox functions (Schwanghart et al., 2014)

% the code below uses 2 types of data that must be obtained from the TopoToolBox functions using the GRIDobj of the topography of the studied area (DEM):
    
    % A STREAMobj including the streams to be studied (computed from DEM): here I stored theses
    % streams in a cell array (eST) which includes 3 streams (eST{1} is stream ST3A, 
    % eST{2} is stream ST3B, and eST{3} is stream ST3C)

    % A is the flow accumulation grid computed from DEM


%% extract relevant data

for i=1:3
    a{i} = getnal(eST{i},A)*DEM.cellsize^2; % extract drained area along each stream
    ghat{i} = 1./(a{i}.^0.45);  % compute the inverse of drained area raised to the power m/n (here taken as 0.45)
    chi{i} = cumtrapz(eST{i},ghat{i}); % compute chi along streams
    z{i}=getnal(eST{i},DEM); % extract the elevation along the streams
    M{i}=horzcat(a{i},ghat{i},chi{i},z{i}); % store the 4 latter parameters in matrix M
end

% make a chi-sorted copy of M

M2{1}=sortrows(M{1},3); % sort matrix of first stream as a function of chi
M2{2}=sortrows(M{2},3); % sort matrix of second stream as a function of chi
M2{3}=sortrows(M{3},3); % sort matrix of third stream as a function of chi


%% Identify confluence

% ST3A-B confluence

for i=1:min([length(M2{1}) length(M2{2})])
    test(i)=M2{1}(i,4)-M2{2}(i,4);
end

IND_CONF=find(test,1,'first'); % extract indice of confluence (in 2nd stream reference)
chi_thresh=M2{1}(IND_CONF,3); % extract chi-value at the confluence

% ST1-ST3 confluence

clear test

for i=1:min([length(M2{1}) length(M2{3})])
    test(i)=M2{1}(i,4)-M2{3}(i,4);
end

IND_CONF_2=find(test,1,'first'); % extract indice of confluence (in 2nd stream reference)
chi_thresh_2=M2{1}(IND_CONF_2,3); % extract chi-value at the confluence


%% compute chi for a range of lost areas scenario

X=0:15000000:300000000; % define grid of total area losses to be tested ST3A+ST3B;
XX=0:1:0; % define grid of area loss to be tested for tributary ST3C;

p(1,:)=[0.5:0.01:1]; % define grid of proportion of lost drained area of ST3A relative to lost drained area of ST3A+ST3B
p(2,:)=1-p(1,:); % proportion of lost drained area of stream ST3B relative to lost drained area of ST3A+ST3B


% compute corrected chi of ST3A and ST3B depending on lost drained areas scenario

for r=1:length(XX)
for i=1:length(X)
    for j=1:length(p)
            for k=1:2 % indice k refer to ST3A or ST3B
                corr_a{r,i,j,k}=M{k}(:,1)+X(i)*(p(1,j)+p(2,j)); % compute drained area along stream k (between the confluences ST3A-B and ST3C) as a function of total lost area (i) and relative proportion of lost area (j)
                corr_a{r,i,j,k}(M{k}(:,3)>chi_thresh)=M{k}(M{k}(:,3)>chi_thresh,1)+X(i)*p(k,j); % update above upper confluence
                corr_a{r,i,j,k}(M{k}(:,3)<=chi_thresh_2)=M{k}(M{k}(:,3)<=chi_thresh_2,1)+X(i)+XX(r); % update below lower confluence
                corr_ghat{r,i,j,k} = 1./(corr_a{r,i,j,k}.^0.45);
                corr_chi{r,i,j,k} = cumtrapz(eST{k},corr_ghat{r,i,j,k});
            end
    end
end
end


%% Project orthogonal to fault

% compute distance from fault (here the trend of the fault is defined as a linear equation in cartesian coordinates: y=mm*x+pp)

mm=-1.4023; % update if needed 
pp=3477817; % update if needed

D{1}=(sqrt((mm*eST{1}.x-eST{1}.y+pp).^2))./(sqrt(mm^2+1)); % compute distance from fault along ST3A
D{2}=(sqrt((mm*eST{2}.x-eST{2}.y+pp).^2))./(sqrt(mm^2+1)); % compute distance from fault along ST3B
D{3}=(sqrt((mm*eST{3}.x-eST{3}.y+pp).^2))./(sqrt(mm^2+1)); % compute distance from fault along ST3C  

D_CONF_12=(sqrt((mm*eST{1}.x(M{1}(:,3)==chi_thresh)-eST{1}.y(M{1}(:,3)==chi_thresh)+pp).^2))./(sqrt(mm^2+1)); % compute distance from fault of the ST3A-B confluence
D_CONF_13=(sqrt((mm*eST{1}.x(M{1}(:,3)==chi_thresh_2)-eST{1}.y(M{1}(:,3)==chi_thresh_2)+pp).^2))./(sqrt(mm^2+1)); % compute distance from fault of the ST3C confluence

% interpolate elevation and corrected chi as a function of fault distance across a regular grid

G=0:2:21000; % define distance from fault as a regular grid (in meters)

[~,G_CONF_IND_12]=min(abs(G-D_CONF_12)); % extract indice of maximum fault distance where ST3A an ST3B do not overlap
G_CONF_IND_12=G_CONF_IND_12-1;

[~,G_CONF_IND_13]=min(abs(G-D_CONF_13)); % extract indice of maximum fault distance where ST3A an ST3C do not overlap
G_CONF_IND_13=G_CONF_IND_13-1;

for r=1:length(XX)
for i=1:length(X)
    for j=1:length(p)
        for k=1:3
            interp_chi{r,i,j,k}=interp1(D{k},corr_chi{r,i,j,k},G); % interpolate corrected chi along streams as a function of distance from fault
        end
    end
end
end
            
for k=1:3
    interp_z{k}=interp1(D{k},M{k}(:,4),G); % interpolate elevation along streams as a function of distance from fault 
end


%% Compute pre-deformation steepness index to make the uplift profile of ST3A and ST3B consistent

dZ1=interp_z{1}(1:G_CONF_IND_12)-interp_z{1}(G_CONF_IND_12); % compute elevation difference of ST3A relative to the ST3A-B confluence
dZ2=interp_z{2}(1:G_CONF_IND_12)-interp_z{2}(G_CONF_IND_12); % compute elevation difference of ST3B relative to the ST3A-B confluence

for r=1:length(XX)
for i=1:length(X)
    for j=1:length(p)
        dChi1{r,i,j}=interp_chi{r,i,j,1}(1:G_CONF_IND_12)-interp_chi{r,i,j,1}(G_CONF_IND_12); % compute chi difference of ST3A relative to the ST3A-B confluence
        dChi2{r,i,j}=interp_chi{r,i,j,2}(1:G_CONF_IND_12)-interp_chi{r,i,j,2}(G_CONF_IND_12); % compute chi difference of ST3B relative to the ST3A-B confluence
    end
end
end

for r=1:length(XX)
for i=1:length(X)
    for j=1:length(p)
        clear mdl
        mdl=fitlm(dChi1{r,i,j}-dChi2{r,i,j},dZ1-dZ2); % compute linear regression between dZ inter-stream difference and dChi inter-stream difference
        SLOPE(r,i,j)=mdl.Coefficients.Estimate(2); % extract slope of the regression (gives the pre-deformation steepness index)
        P_VAL(r,i,j)=mdl.Coefficients.pValue(2); % extract p-value of the regression
        R_SQ(r,i,j)=mdl.Rsquared.Ordinary; % extract R-square of the regression
    end
end
end


SLOPE(P_VAL>0.05)=NaN; % remove scenarios that have p-values above 0.05
SLOPE(R_SQ<0.5)=NaN; % remove scenarios that have R-square below 0.5


%% compute slope of chi-plot between the 2 confluences (main trunk T)

for r=1:length(XX)
for i=1:length(X)
    for j=1:length(p)
                clear mdl
                mdl=fitlm(interp_chi{r,i,j,1}(G_CONF_IND_12:G_CONF_IND_13),interp_z{1}(G_CONF_IND_12:G_CONF_IND_13)); % compute linear model for the slope of the chi-plot between the confluence and 6000 m from the fault
                MDL(r,i,j)=mdl.Coefficients.Estimate(2); % store slope of chi-plot between the 2 confluences 
                P_VALUE(r,i,j)=mdl.Coefficients.pValue(2); % store p value of the regression
                MISFIT(r,i,j)=mdl.Coefficients.Estimate(2)/SLOPE(r,i,j); % compute ratio between slope of chi-plot yielded by co-tectonic lines and that obtained between the 2 confluences
    end
end
end

MISFITc=MISFIT; % make a copy of MISFIT matrix

% remove solutions which do not fit with independent range of pre-deformation chi-plot

for r=1:length(XX)
for i=1:length(X)
    for j=1:length(p)
            if SLOPE(r,i,j)>=BSmed-BSuncert & SLOPE(r,i,j)<=BSmed+BSuncert % here BSmed and BS uncert refer to independent constraints on the pre-beheading steepness index: 29.2+/-7.0  (BSmed=29.2; BSuncert=7.0)
                MDL(r,i,j)=MDL(r,i,j);
                P_VALUE(r,i,j)=P_VALUE(r,i,j);
                MISFIT(r,i,j)=MISFIT(r,i,j);
            else
                MDL(r,i,j)=NaN;
                P_VALUE(r,i,j)=NaN;
                MISFIT(r,i,j)=NaN;
            end
    end
end
end


%% Create logical matric that contains NaN for all non-viable solutions

LOGIC=MISFIT;

for r=1:length(XX)
for i=1:length(X)
    for j=1:length(p)
        if MISFIT(r,i,j)<1.0
            LOGIC(r,i,j)=NaN;
        end
    end
end
end


%% Compute uplift


for r=1:length(XX)
for i=1:length(X)
    for j=1:length(p)
        if isnan(LOGIC(r,i,j))==0
            clear kk1
            [~,kk1]=max(interp_chi{r,i,j,1}); % extract indice of the upper termination of ST3A
            DELTA_ZZ{r,i,j}=SLOPE(r,i,j).*interp_chi{r,i,j,1}+interp_z{1}(kk1)-SLOPE(r,i,j).*max(interp_chi{r,i,j,1})-interp_z{1}; % extract uplift difference along ST3A relative to its upper termination
            REL_UP{r,i,j}=-DELTA_ZZ{r,i,j}+DELTA_ZZ{r,i,j}(G_CONF_IND_13); % compute uplift relative to the ST3C confluence
        else
            DELTA_ZZ{r,i,j}=NaN;
            REL_UP{r,i,j}=NaN;
        end
    end
end
end



%% Make contour plot for inversion (Fig. 6 of the article)

vvvvv=[1 1];
vvv=[1.1 1.2 1.3 1.4 1.5];
vv=[BSmed-BSuncert BSmed+BSuncert];
v=[BSmed BSmed];
MM=contour(reshape(SLOPE,[length(X),length(p)]),vv);

level1=MM(:,2:MM(2,1)+1);
level2=MM(:,MM(2,1)+3:end);

x=[level1(1,:) flip(level2(1,:))];
y=[level1(2,:) flip(level2(2,:))];

figure;
fill(x,y,[0.8 0.8 0.8])
hold on
contour(reshape(SLOPE,[length(X),length(p)]),v,'Color','k','LineWidth',3)
contour(reshape(MISFITc,[length(X),length(p)]),vvvvv,'Color','b','LineWidth',3)
contour(reshape(MISFITc,[length(X),length(p)]),vvv,'Color','b','LineWidth',0.5,'LineStyle','-.')

axis([1 51 2 21])

scat=scatter(35.12,12.333,'*')
scat.LineWidth=2
scat.SizeData=100
scat.MarkerEdgeColor='r'

xlabel('ST3A catchment proportion','FontSize',18,'FontWeight','bold','Interpreter','Latex')
ylabel('Total area lost ($km^2$)','FontSize',18,'FontWeight','bold','Interpreter','Latex')

xticks([1 11 21 31 41 51]);
xticklabels({'0.5','0.6','0.7','0.8','0.9','1.0'})
yticks([1 7 13 19]);
yticklabels({'0','90','180','270'})

plot([-1 -2],[-1 -2],'Color','b','LineWidth',1)
plot([-1 -2],[-1 -2],'Color','b','LineWidth',0.5,'LineStyle','-.')

legend('Steepness index = 29.2±7.0','',['Steepness index' newline 'misfit = 0%'],'','','',['Steepness index' newline 'misfit = 10 to 50%'],'','FontSize',11)

ta=annotation('textarrow',[0.7755 0.67],[0.30 0.53],'String',{'observed','value'}','Interpreter','Latex')
ta.Color='r'
ta.FontSize=18

exportgraphics(gca,'NewFig6.jpg','resolution',1200)



%% Plot uplift (Fig.7 and S4 of the article)

% Prepare data

Temp1=[nan];
Temp2=[nan];
Temp3=[nan];
Temp4=[nan];

for r=1:length(XX)
for i=1:length(X)
    for j=1:length(p)
        if sum(sum(sum(isnan(REL_UP{r,i,j}))))<numel(REL_UP{r,i,j})
            Temp1=horzcat(Temp1,REL_UP{r,i,j});
            Temp2=horzcat(Temp2,G);
            Temp3=horzcat(Temp3,repmat(p(1,j),1,numel(G)));
            Temp4=horzcat(Temp4,repmat(X(1,i),1,numel(G)));
        end
    end
end
end

Temp2(isnan(Temp1)==1)=[];
Temp3(isnan(Temp1)==1)=[];
Temp4(isnan(Temp1)==1)=[];
Temp1(isnan(Temp1)==1)=[];
TempT=vertcat(Temp1,Temp2,Temp3,Temp4);

cols=size(TempT,2);
P=randperm(cols);
TempT=TempT(:,P);

BestTempT=TempT(:,TempT(3,:)==0.85);
BestTempT=BestTempT(:,BestTempT(4,:)==165000000);

% Plot color-coded with ST3A proportion

figure;
plot([0 G(G_CONF_IND_13)],[0 0],'k','LineStyle','-.','LineWidth',1.5)
hold on
scatter(TempT(2,2:end),TempT(1,2:end),1,TempT(3,2:end))
colorbar
axis([0 G(G_CONF_IND_13) -10 50])
plot([G(G_CONF_IND_12) G(G_CONF_IND_12)],[-10 50],'b','LineStyle','-.','LineWidth',0.75)

xlabel('Distance from fault (m)','FontSize',18,'FontWeight','bold','Interpreter','Latex')
ylabel('Relative uplit (m)','FontSize',18,'FontWeight','bold','Interpreter','Latex')

xticks([0 1000 2000 3000 4000]);
xticklabels({'0','1000','2000','3000','4000'})
yticks([0 10 20 30 40 50]);
yticklabels({'0','10','20','30','40','50'})

cb1=colorbar;
set(cb1,'Position',[0.32 0.30 0.02 0.2])
ylabel(cb1,{'ST3A catchment','proportion'},'FontSize',14,'Interpreter','Latex')
cb1.Label.Position=[-2.8 0.865 0]

ta=annotation('textarrow',[0.6 0.7],[0.82 0.82],'String',{'ST3A-B   ','confluence   '}','Interpreter','Latex')
ta.Color='b'
ta.FontSize=14

exportgraphics(gca,'NewFig7.jpg','resolution',1200)


% Plot color-coded with total loss area

figure;
plot([0 G(G_CONF_IND_13)],[0 0],'k','LineStyle','-.','LineWidth',1.5)
hold on
scatter(TempT(2,2:end),TempT(1,2:end),1,TempT(4,2:end))
colorbar
axis([0 G(G_CONF_IND_13) -10 50])
plot([G(G_CONF_IND_12) G(G_CONF_IND_12)],[-10 50],'k','LineStyle','-.','LineWidth',0.75)

xlabel('Distance from fault (m)','FontSize',18,'FontWeight','bold','Interpreter','Latex')
ylabel('Relative uplit (m)','FontSize',18,'FontWeight','bold','Interpreter','Latex')

xticks([0 1000 2000 3000 4000]);
xticklabels({'0','1000','2000','3000','4000'})
yticks([0 10 20 30 40 50]);
yticklabels({'0','10','20','30','40','50'})

cb1=colorbar;
set(cb1,'Position',[0.32 0.34 0.02 0.2])
ylabel(cb1,'Total area loss ($km^2$)','FontSize',14,'Interpreter','Latex')
cb1.XTick=[150000000 200000000 250000000 300000000]
cb1.XTickLabel={'150','200','250','300'}



exportgraphics(gca,'FigSupp4.jpg','resolution',1200)



