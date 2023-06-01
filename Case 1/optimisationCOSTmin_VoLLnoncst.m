function [ENS,LOL,load,P_value,P_stor_value,E_stor_value,year_MC] = optimisationCOSTmin_VoLLnoncst( ...
    ccgt_available_power,wind_available_power,TS_load,VOLL,depth_voll,columns,stor_cap,stor_ener,Nt,eta_stor,cost)
% This function returns different arrays of the results of the adequacy
% assessment, such as the hourly energy not served or the energy level of
% the storage when the objective function is to minimise total cost, with a
% non constant value of lost load.
%
% Input:
%   - ccgt_available_power: vector of available power of classical units (in this
%   case, ccgt units) [MW]
%   - wind_available_power: matrix of hourly available RES power of (in this
%   case, offshore wind) for 38 climate years (1982 to 2019) [MW]
%   - TS_load: matrix containing hourly loads for 35 climate years (1982 to
%   2016) [MW]
%   - VOLL: value of lost load, corresponding to the cost of one MWh of
%   energy not supplied [€/MWh]
%   - depth_voll: depth of the load shedding for which a constant VoLL is
%   assumed [MW]
%   - columns: number of columns to store the load shedding [-]
%   - stor_cap: capacity of storage [MW]
%   - stor_ener: size of storage [MWh]
%   - Nt: time length of the adequacy assessment [-]
%   - eta_stor: one-way efficiency of storage [-]
%   - cost: matrix containing cost of different production technologies
%   considered [€/MWh]
%
% Output:
%   - ENS: 8760x1 vector equal to 0 when there is no shortage and equal to
%   the value of the storage [MW]
%   - LOL: 8760x1 vector equal to 0 when ENS<0.5 and 1 when ENS>=0.5 [-]
%   - load: 8760x1 vector of the load for the considered year [MW]
%   - P_value: 8760x2 matrix containing non RES and RES power available for
%   the considered year [MW]
%   - P_stor_value: 8760x2 matrix containing in the 1st column the discharging 
%   power while in the 2nd column the charging power of storage
%   - E_stor_value: 8760x1 vector containing the energy level of storage
yalmip('clear');
ones_Nt = ones(1,Nt);
VOLL_vector = linspace(VOLL,2*VOLL,columns)'; %every step, there is a increase of the VOLL

year_MC = randi(35,1); %generates random number to select already-made TS for load
Pmax = [ccgt_available_power wind_available_power(:,year_MC)];
load = TS_load(:,year_MC);

P = sdpvar(Nt,2,'full');
P_stor_prod = sdpvar(Nt,1,'full');
P_stor_cons = sdpvar(Nt,1,'full');
E_stor = sdpvar(Nt,1,'full');
% This LS_matrix will allow to 'sort' the load shedding according to its
% depth e.g., a 12 GW load shedding will cost more than 50% of a 8 GW one
% (if VoLL was constant, cost of 12 GW LS = 1.5*(8 GW LS) (LS stands for
% load shedding)
LS = sdpvar(Nt,1,'full');
LS_matrix = sdpvar(Nt,columns,'full'); 
y = sdpvar(1,columns,'full');

%constraints: power should be lower or equal to available power and >= 0
Constraints = [P <= Pmax, P >= zeros(Nt,2)];

%Storage constraints
%   Max and minimum power constraints; max energy constraints
Constraints = [Constraints, P_stor_cons >= 0; P_stor_cons <= stor_cap];
Constraints = [Constraints, P_stor_prod >= 0; P_stor_prod <= stor_cap];
Constraints = [Constraints, E_stor >= 0; E_stor <= stor_ener];
%   initial level
Constraints = [Constraints, E_stor(1,1) == stor_ener/2];
%   efficiency of storage
Constraints = [Constraints, E_stor(2:end,1) == E_stor(1:end-1,1) - P_stor_prod(1:end-1,1)*(1/eta_stor) + P_stor_cons(1:end-1,1)*eta_stor];

yalmip('nvars')
%LS = max(load-sum(P,2)-P_stor_prod+P_stor_cons,0);
Constraints = [Constraints, LS>=load-sum(P,2)-P_stor_prod+P_stor_cons,LS>=0];
% F = implies(load-sum(P,2)-P_stor_prod+P_stor_cons<=0,LS==0);
% Constraints = [Constraints, F, LS <= 100e3, LS>=0];


% Ensures that the rows of LS_matrix corresonds to the hourly load shedding
Constraints = [Constraints, sum(LS_matrix,2)==LS];
% Each element of the matrix must be between 0 and depth_voll MW (depth_voll 
% being the arbitrary chosen step, e.g. a 2.1 GW load shedding would be 
% split into 3 parts (1+1+0.1)
Constraints = [Constraints; LS_matrix >=0, LS_matrix <= depth_voll];
Constraints = [Constraints, y == ones_Nt*LS_matrix];

% Objective function (minimisation of cost)
tic
Objective = ones_Nt*P*cost' + y*VOLL_vector;
toc

% Set some options for YALMIP and solver
options = sdpsettings('verbose',0,'solver','cplex','cplex.timelimit',120, 'cplex.mip.tolerances.mipgap',1e-4);
sol = optimize(Constraints,Objective,options)

P_value(:,1) = value(P(:,1)); P_value(:,2) = value(P(:,2));
P_stor_value(:,1) = value(P_stor_prod); P_stor_value(:,2) = value(P_stor_cons);
E_stor_value = value(E_stor);
PNS = value(LS);
ENS = max(PNS,0);
LOL = round(min(ENS,1)); %to avoid computational errors
end