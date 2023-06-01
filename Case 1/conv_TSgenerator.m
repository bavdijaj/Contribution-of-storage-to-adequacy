function [conv_TS] = conv_TSgenerator(FO_una,FO_MTTR,PO_days)
% This function gives yearly availability of a unit according to its forced
% and planned outage (resp. FO & PO) rates and mean time to repair (MTTR)
%
% Input:
%   - FO_una: Annual availability due to forced outages, in %
%   - FO_MMTR: Mean time to repair after a forced outage, in days
%   - PO_days: number of days during which planned outages occur
% Output:
%   - conv_TS: returns a 8760x1 vector with values between 0 and 1 which
%   corresponds to hourly availabilities
conv_TS = [];
avai_inf = 1-FO_una; %asymptotic availability corresponding to the % of the time the unit is available
FO_MTTF = avai_inf*FO_MTTR*24/(1-avai_inf);
FO_fail_rate = 1/FO_MTTF; FO_rep_rate = 1/(FO_MTTR*24);
AU1 = 1; %suppose all units are working 

while length(conv_TS) < 8760
%    My Technique
% fail = binornd(1,FO_fail_rate);
% if fail == 0
%     conv_TS(end+1,1) = 1;
%     %i = i+1;
% else
%     conv_TS(end+1:end+FO_MTTR*24,1) = 0;
%     %i = i + 24;
% end

%Chiara's technique
U1= rand;
if AU1 == 1 %si l unit1 est disponible
    lambda1 = FO_fail_rate;
    T1 = round(-1/(lambda1)*log(U1)); % Formule pour calculer le temps pdt lequel Unit1 est disponible
    conv_TS(end+1:end+T1) = 1; %la valeur de la fonction pendant cette durée vaut 1 (Unit1 disponible)
    AU1 = 0; %une fois que cette période est passé, l'unit1 n'est plus disponible
else
    lambda1 = FO_rep_rate;
    T1 = round(-1/(lambda1)*log(U1));
    conv_TS(end+1:end+T1) = 0;
    AU1 = 1;
end

end

%just to allow the placement of PO to happen without error; it will have no
%impact on the result as it is troncated at the end of the function
conv_TS(end+1:end+200) = 1; 


margin = 2;
day = 1.5;
time_step = round(365/day)-margin;
po_moments_days = randperm(time_step,PO_days); %will spread POs randomly with at least day days between them
po_moments = po_moments_days*(day*24);
for i=1:PO_days
    j = po_moments(i)+randi([-12 12]); %so as not to have exactly the same profile for every unit
    while isequal(isequal(conv_TS(j:j+23),ones(1,24)),0) %checks if any forced outages had already happened
        j = j+1;
    end
    conv_TS(j:j+24)=0;
end

conv_TS = conv_TS(1:8760); %only keeping the TS for one year, i.e., 8760 hours
end