clear all; close; clc;

Nt = 8760; %number of time steps considered
%% Definition of the strategy
% Use of a variable, named stategy, to define the strategy of the storage
% strategy = 1 --> Cost minimisation, with constant VOLL
% strategy = 2 --> Cost minimisation, with non constant VOLL
% strategy = 3 --> LOLE minimisation, with constant VOLL
strategy = 3;
%% Data import & initialisation
%installed capacity of RES
wind_off = 109.9727e3; %offshore wind capacity in MW,
cost_windOff = 0;

%installed storage
stor_ener = 1*36.74e3; %storage size in MWh
stor_cap = 1*9.5e3; %storage capacity in MW
eta_stor = 0.9; %one-way efficiency of storage --> ASSUMPTION to confirm

%installed conventional (CCGT present 1 in PEMMDB National Estimates 2022)
ccgt_units = 100; %number of units installed
ccgt_size = 37551.28/ccgt_units; %capacity of ccgt unit in MW
cost_ccgt = 49.8; %€/MWh
ccgt_FO_rate = 0.05; % ATTENTION: be sure how to modelise this value
ccgt_FO_MTTR = 1;
ccgt_PO_rate = 27; % number of days during which POs occur
available_ccgt = zeros(Nt,ccgt_units);
 
%historical capacity factor 
TS_windOff = table2array(readtable("2022 National Estimates/2022_UK011_windOffshore_CF.csv"));
wind_available_power = wind_off*TS_windOff;

cost = [cost_ccgt cost_windOff];

% historical load
TS_load = table2array(readtable("2022 National Estimates/2022_UK_load_2030estimated.csv"));

VOLL = 17000; %Value of lost load, in €/MWh 
%% Optimisation 
treshold = 0.001; %same treshold as Elia used
alpha = 1; %initialisation of the coefficient of variation
criteria = 1; %initialisation of the increment coefficient of variation
LOL_yearly = zeros(1000,1); %array in which different LOLEs will be stored
LOL_hourly = zeros(8760,1000);
ENS_yearly = zeros(1000,1);
ENS_hourly = zeros(8760,1000);
Load_hourly = zeros(8760,1000); P_yearly = zeros(8760,1000,2); P_stor_yearly = zeros(8760,1000,2);
E_stor_yearly = zeros(8760,1000);
LS_hourly = zeros(8760,1000);
years = zeros(1000,1);
nonres_available_power = zeros(Nt,1000);

inc_yearly = ones(1000,1);
i = 1;

rng('default'); %to reinitialise the random number generator such as to have reproducible simulations
sim_iter = 200; %number of iterations to perform, same size for each strategy to be able to compare results
                %200 was deemed to be sufficient from results of previous
                %simulations
while i < sim_iter
%while criteria > treshold
    for j=1:ccgt_units
        % for each ccgt, we generate an hourly TS of the available power
        available_ccgt(:,j) = ccgt_size*conv_TSgenerator(ccgt_FO_rate,ccgt_FO_MTTR,ccgt_PO_rate);
    end
    %sums of each row, giving available power for each hour
    ccgt_available_power = sum(available_ccgt,2);
    nonres_available_power(:,i) = ccgt_available_power;
    %plot(ccgt_available_power); hold on;
    if strategy == 1
        [ENS_hourly(:,i), LOL_hourly(:,i), Load_hourly(:,i), P_yearly(:,i,:), P_stor_yearly(:,i,:), E_stor_yearly(:,i), years(i)] = optimisation(ccgt_available_power, ...
            wind_available_power,TS_load,VOLL,stor_cap,stor_ener,Nt,eta_stor,cost);
    elseif strategy == 2
        depth_voll = 1e3;%depth of the load shedding for which a constant VoLL is assumed
        
        %to determine numbers of columns in LS_matrix; 60 because this is
        %around the max pea.7k, so to ensure there is enough columns to store
        %the load shedding
        columns = 60e3/depth_voll; 
        [ENS_hourly(:,i), LOL_hourly(:,i), Load_hourly(:,i), P_yearly(:,i,:), P_stor_yearly(:,i,:), E_stor_yearly(:,i),years(i)] = optimisationCOSTmin_VoLLnoncst(ccgt_available_power, ...
            wind_available_power,TS_load,VOLL,depth_voll,columns,stor_cap,stor_ener,Nt,eta_stor,cost);
    elseif strategy == 3
        [ENS_hourly(:,i), LOL_hourly(:,i), Load_hourly(:,i), P_yearly(:,i,:), P_stor_yearly(:,i,:), E_stor_yearly(:,i), LS_hourly(:,i),years(i)] = optimisationLOLEmin(ccgt_available_power, ...
            wind_available_power,TS_load,stor_cap,stor_ener,Nt,eta_stor);
    end

    ENS_yearly(i) = sum(ENS_hourly(:,i));
    LOL_yearly(i) = sum(LOL_hourly(:,i));
    
    EENS = sum(ENS_yearly)/i;
    var_EENS = var(ENS_yearly(ENS_yearly>0));
    inc = abs(sqrt(var_EENS)/EENS-alpha)/alpha;
    alpha = sqrt(var_EENS)/EENS;
    inc_yearly(i) = inc;
    if i > 4
        criteria = mean(inc_yearly(i-4:i)) %to avoid only one small value stopping 
    else 
        criteria = inc
    end
    i = i+1
    clear year_MC EENS var_EENS; clear conv_TSgenerator optimisation;
    yalmip('clear');
end

[lol_cdf, lol_dist] = ecdf(LOL_yearly(LOL_yearly>0));
LOLE = mean(LOL_yearly(LOL_yearly>0));

%% Test 
Real_LS = max(Load_hourly-P_yearly(:,:,1)-P_yearly(:,:,2)-P_stor_yearly(:,:,1)+P_stor_yearly(:,:,2),0);

%% 
EENS = mean(ENS_yearly(ENS_yearly > 0));
ENS_only = ENS_yearly(ENS_yearly>0);
LOL_only = LOL_yearly(LOL_yearly>0);
std_LOLE = std(LOL_only)/sqrt((sim_iter-1));
std_EENS = std(ENS_only)/sqrt((sim_iter-1));

%% visualisation and post-processing
figure(1)
MC_year1 = 11; %select which MC year you want to show
sample = 2;
tiledlayout(3,1,'Padding','Compact')
nexttile
plot(P_yearly(end-24*sample:end,MC_year1,1)/1e3,LineWidth=1.2);
hold on;
plot(P_yearly(end-24*sample:end,MC_year1,2)/1e3,LineWidth=1.2);
plot(Load_hourly(end-24*sample:end,MC_year1)/1e3,LineWidth=1.2);
title('Generation and demand curves'); ylabel('Power [GW]'); xlabel('Time [h]');
legend('CCGT','Offshore wind','Load','FontSize',10.5)
nexttile

yyaxis left
plot(E_stor_yearly(end-24*sample:end,MC_year1)/1e3,LineWidth=1.2);
ylabel('Energy level [GWh]')
hold on;

yyaxis right
h = zeros(1,2);
h(1) = plot(P_stor_yearly(end-24*sample:end,MC_year1,1)/1e3,LineWidth=1.2, DisplayName='Discharging');
h(2) = plot(P_stor_yearly(end-24*sample:end,MC_year1,2)/1e3,LineWidth=1.2, DisplayName='Charging');
title('Storage curves');
xlabel('Time [h]'); ylabel('Power [GW]');
legend(h,'FontSize',10.5);

nexttile
yyaxis left
stem(LOL_hourly(end-24*sample:end,MC_year1),'filled',LineStyle=':');
yticks([0 1])
ylabel('LOL [-]')
hold on;

yyaxis right
plot(ENS_hourly(end-24*sample:end,MC_year1)/1e3,LineWidth=1.2)
title('Loss of Load'); ylabel('Load shedding [GW]'); xlabel('Time [h]');

%%
figure(2)
MC_year2 = 13; %select which MC year you want to show
tiledlayout(3,1)
nexttile
plot(P_yearly(end-24*4:end,MC_year2,1)/1e3,LineWidth=1.2);
hold on;
plot(P_yearly(end-24*4:end,MC_year2,2)/1e3,LineWidth=1.2);
plot(Load_hourly(end-24*4:end,MC_year2)/1e3,LineWidth=1.2);
title('Generation and demand curves'); ylabel('Power [GW]'); xlabel('Time [h]');
legend('CCGT','Offshore wind','Load')
nexttile

plot(E_stor_yearly(end-24*4:end,MC_year2)/1e3,LineWidth=1.2);
hold on;
plot(P_stor_yearly(end-24*4:end,MC_year2,1)/1e3,LineWidth=1.2);
plot(P_stor_yearly(end-24*4:end,MC_year2,2)/1e3,LineWidth=1.2);
title('Storage curves');
xlabel('Time [h]'); ylabel('Power [GW]');
legend('Energy level [GWh]', 'Discharging', 'Charging');

nexttile
plot(LOL_hourly(end-24*4:end,MC_year2));
title('Loss of Load'); ylabel('LOL [-]'); xlabel('Time [h]');

%% 
figure(3)

histogram(lol_dist,'BinWidth',25,'Normalization','probability');
%title("Probability distribution of LOLE (cost min. with non constant VoLL)")
xlabel("LOLE [h]"); ylabel("Probability [-]")

text_avg = append(num2str(round(LOLE)),' hours');
xline(LOLE,'-',{'Average',text_avg},'LabelHorizontalAlignment','left');
%%
figure(4)
ecdf(lol_dist); 
%title('Cumulative distribution of LOLE (cost min. with non constant VoLL)')
xlabel("LOLE [h]");ylabel('Percentile [-]');
xline(LOLE,'-',{'Average',text_avg},'LabelHorizontalAlignment','left');
%% 
i=1;
while i < sim_iter
%while criteria > treshold
    for j=1:ccgt_units
        % for each ccgt, we generate an hourly TS of the available power
        available_ccgt(:,j) = ccgt_size*conv_TSgenerator(ccgt_FO_rate,ccgt_FO_MTTR,ccgt_PO_rate);
    end
    %sums of each row, giving available power for each hour
    ccgt_available_power = sum(available_ccgt,2);
    nonres_available_power(:,i) = ccgt_available_power;
    i=i+1;
end
avg = mean(mean(nonres_available_power))
stddev = std(mean(nonres_available_power,2))