function Out = MainFunction(In)

% Read xls Data
sysData = xlsread('ieee118cdf-ee252.xls'); 
busData = csvread('Bus_Data.csv',1,0);
error = 0.001;                                %% Acceptable error in kW or kVar

% %%%%%%%%%%%%%%% Build YBus Using Building Block Method %%%%%%%%%%%%%%%%%
yBus = BuildingBlockMethod(sysData);

fName = 'yBus.csv';      
dlmwrite(fName,yBus,'delimiter',',','precision','%0.5f');  

Vbus = runFDPF_XB(busData,yBus,sysData,error);
%Vbus(:,3) = Vbus(:,3).*(180/pi);               %% TO DELETE: Convert angle to degrees

fName = 'Bus_Voltage.csv';      
dlmwrite(fName,Vbus,'delimiter',',','precision','%0.5f');  

end

function bBlock_yBus = BuildingBlockMethod(SystemData)
[numBus,lRow] = GetNumberOfBusesRows(SystemData);
yBus = zeros(numBus);
bBlock = [1 -1; -1 1];              %Building Block

for i = 1:lRow
    fromBus = SystemData(i,1);
    toBus = SystemData(i,2);
    ySeries = 1/(SystemData(i,4) + (1i) * SystemData(i,5));
    yShunt = 0.5*(1i) *SystemData(i,6);
    
    %Add Series Admittance
    subMatrix = bBlock .* ySeries;
    yBus(fromBus,fromBus) = yBus(fromBus,fromBus) + subMatrix(1,1);
    yBus(fromBus,toBus) = yBus(fromBus,toBus) + subMatrix(1,2);
    yBus(toBus,fromBus) = yBus(toBus,fromBus) + subMatrix(2,1);
    yBus(toBus,toBus) = yBus(toBus,toBus) + subMatrix(2,2);
    
    %Add shunt Admittance
    subMatrix = 1 * yShunt;
    yBus(fromBus,fromBus) = yBus(fromBus,fromBus) + subMatrix(1,1);
    
    %Add shunt Admittance
    subMatrix = 1 * yShunt;
    yBus(toBus,toBus) = yBus(toBus,toBus) + subMatrix(1,1);
    
end

bBlock_yBus = yBus;

end

function [numBus,lRow] = GetNumberOfBusesRows(systemData)

lRow = size(systemData,1);
numBus = 0;

for i = 1:lRow
   for y = 1:2
       if systemData(i,y) > numBus
          numBus = systemData(i,y);
       end
   end
end

end


function Vbus = runFDPF_XB(busData,yBus,branchData,error)

nBus = size(busData,1);                                             % Number of buses in the System
Vbus = zeros(nBus,3);                                               % Bus voltages

[PQbus,PVbus,SLACKbus] = GroupBuses(busData);                       % Group Buses based on Bus Type
[B1,B2] = getB1B2(yBus,busData,branchData,PQbus,PVbus,SLACKbus);    % Get B', B"

fName = 'B1.csv';      
dlmwrite(fName,B1,'delimiter',',','precision','%0.5f');
fName = 'B2.csv';      
dlmwrite(fName,B2,'delimiter',',','precision','%0.5f');

% Known P's
nPsch = busData((busData(:,2) < 3),1);                              % Buses w/ known Ps
Psch = (busData(nPsch,7) - busData(nPsch,5))/100;                   % Injected Power = Gen Power - Load Power, w/ 100MVA base

% Known Q's
Qsch = (busData(PQbus,8) - busData(PQbus,6))/100;                   % Injected Reactive Power = Gen Reactive Power - Load Reactive Power w/ 100MVA base

%% Initialize Calculated Qs and Ps
Qcalc = zeros(size(PQbus,1),1);
Pcalc = zeros(size(nPsch,1),1);

%% Determine/Initialize Values of the Bus Voltages
Vbus(:,1) = busData(:,1);
% bus type is slack bus
Vbus(SLACKbus,2) =  busData(SLACKbus,3);
Vbus(SLACKbus,3) = (busData(SLACKbus,4) * pi)/ 180;     %% convert angle to radians         
% bus type is PV bus
Vbus(PVbus,2) = busData(PVbus,3);                       
Vbus(PVbus,3) = (busData(SLACKbus,4) * pi)/ 180;        %% Initialize angle to slack bus angle
% bus type is PQ bus
Vbus(PQbus,2) = busData(SLACKbus,3);                    %% Initialize voltage to voltage of slack bus 
Vbus(PQbus,3) = (busData(SLACKbus,4) * pi)/ 180;        %% Initialize angle to slack bus angle

MaxError = 1000000000000000000;                         %% Initialize Maximum error to a large value
iteration = 0;
%% Perform Load Flow
while MaxError > error
    %% Calculate P
    %% Pi = sum( |Yij||Vi||Vj|cos(Aij + Aj - Ai) )
    for i = 1:size(nPsch,1)
        Pcalc(i,1) = 0;
        for j = 1:nBus
            r = nPsch(i,1);     % Corresponding Bus Number
            c = j;              % Corresponding Bus Number
            
            Yij = abs(yBus(r,c));   
            Aij = angle(yBus(r,c));  
            Vi = Vbus(r,2);
            Vj = Vbus(c,2);
            Ai = Vbus(r,3);
            Aj = Vbus(c,3);
            
            Pcalc(i,1) = Pcalc(i,1) + Yij*Vi*Vj*cos(Aij+Aj-Ai);
        end
    end

    %% Calculate Q
    %% Qi = sum( -|Yij||Vi||Vj|sin(Aij + Aj - Ai) )
    for i = 1:size(PQbus,1)
        Qcalc(i,1) = 0;
        for j = 1:nBus
            r = PQbus(i,1);     % Corresponding Bus Number
            c = j;              % Corresponding Bus Number
            
            Yij = abs(yBus(r,c));   
            Aij = angle(yBus(r,c));  
            Vi = Vbus(r,2);
            Vj = Vbus(c,2);
            Ai = Vbus(r,3);
            Aj = Vbus(c,3);
            
            Qcalc(i,1) = Qcalc(i,1) - Yij*Vi*Vj*sin(Aij+Aj-Ai);
        end
    end
    
   
    %% Update Values of Voltage Magnitudes and angles
    
    dP = (Psch - Pcalc)./Vbus(nPsch,2);    % dP = [Pactual - Pcalculated(k)] / |V| 
    dAngle = B1\dP;                        % dP/|V| = B' * dAngle or dAngle = inv(B') * (dP/|V|)
    
    dQ = (Qsch - Qcalc)./Vbus(PQbus,2);    % dQ = [Qactual - Qcalculated(k)] / |V| 
    dV = B2\dQ;                            % dQ/|V| = B" * d|V|
    
    Vbus(nPsch,3) = Vbus(nPsch,3) + dAngle;
    Vbus(PQbus,2) = Vbus(PQbus,2) + dV;
    
    %% Compute for Error
    Perror = (Psch-Pcalc)*100*1000; %in kW
    Qerror = (Qsch-Qcalc)*100*1000; %in kVar
    
    MaxError = 0;
    for i = 1:size(Perror,1)
        if abs(Perror(i,1)) > MaxError
            MaxError = abs(Perror(i,1));
        end
    end
    for i = 1:size(Qerror,1)
        if abs(Qerror(i,1)) > MaxError
            MaxError = abs(Qerror(i,1));
        end
    end
 
    iteration  = iteration +1;
       
end
    fprintf('Converged after %i iterations with an error of %f kW/kVar.\n',iteration,MaxError); 
%     fprintf('\nBus Number |MW Actual| MW Calculated|kW error\n')
%     display([nPsch,Psch*100,Pcalc*100,(Psch-Pcalc)*100*1000]);
%     fprintf('\nBus Number |MVar Actual| MVar Calculated|kVar error\n')
%     display([PQbus,Qsch*100,Qcalc*100,(Qsch-Qcalc)*100*1000]);
    
end

function [B1,B2] = getB1B2(yBus,busData,branchData,PQbus,PVbus,SLACKbus)
  
    
    B1 = zeros(size(yBus,1)-1);                     % size of B' = number of Buses - Slack bus
    
    B2 = zeros(size(PQbus,1));                      % size of B" = number of PQ Buses
    
    nPsch =  busData((busData(:,2) < 3),1);         % Buses w/ known Ps
    
    %% Compute elements of B'
    %% Fill the off-diagonal elements of B1
    for i = 1:size(B1,1)
        for j = 1:size(B1,1)
            if j ~= i 
                r = nPsch(i,1);         % corresponding Bus Number
                c = nPsch(j,1);         % corresponding Bus Number
                
                Bij = imag(-yBus(r,c)); % Imaginary Component of the Negative element of Ybus
                B1(i,j) = Bij;          % Bij' = Bij
            end
        end
    end
    
    %% Fill the diagonal elements B1
    for i = 1:size(B1,1)
        for j = 1:size(B1,1)
            if j ~=i
                B1(i,i) = B1(i,i) - B1(i,j);    % Bii' = -sum(Bij') except for slack Bus
            end
        end
        % Determine Bij' for Slack Bus
        r = nPsch(i,1);                         % corresponding Bus Number
        c = SLACKbus;                           % corresponding Bus Number
        Bij = imag(-yBus(r,c));                 % Imaginary Component of the Negative element of Ybus
        
        B1(i,i) = B1(i,i) - Bij;                % Bii' = -sum(Bij')
        
    end
    
    %% Compute Elements of B" %%%%%%%%
    %% Fill the off-diagonal elements of B"
    for i = 1:size(B2,1)
        for j = 1:size(B2,1)
            if j ~= i 
                r = PQbus(i,1);             % corresponding Bus Number
                c = PQbus(j,1);             % corresponding Bus Number
                
                Bij = imag(0-yBus(r,c));     % Imaginary Component of the Negative element of Ybus
                Gij = real(0-yBus(r,c));     % Real Component of the Negative element of Ybus
                
                if Bij ~= 0
                    B2(i,j) = ((Bij^2) + (Gij^2)) /Bij ;          % Bij" = (Bij ^2 + Gij^2) / Bij
                end
            end
        end
    end
    
%% Fill the diagonal elements B2
    for i = 1:size(B2,1)
        
        for j = 1:size(yBus,1)
            r = PQbus(i,1);                             % corresponding Bus Number
            c = j;                                      % corresponding Bus Number
            if r ~=c
                Bij = imag(-yBus(r,c));                 % Imaginary Component of the Negative element of Ybus
                Gij = real(-yBus(r,c));                 % Real Component of the Negative element of Ybus
                B2ij = ((Bij^2) + (Gij^2)) /Bij ;

                if Bij ~= 0 
                    B2(i,i) = B2(i,i) - B2ij ;          % Bii" = -sum(Bij")
                end
            end
        end
        % Determine Value of Equivalent shunt Reactance
        Bi0 = GetBi0(branchData,r);
        
        B2(i,i) = B2(i,i) -2*Bi0;                       % Bii" = -sum(Bij") -2Bi0
        
    end
    
end

%% Get the shunt Reactance
function Bi0 = GetBi0(branchData,busNumber)
    Bi0 = 0;
    
    for i = size(branchData,1)
        if (branchData(i,1) == busNumber) || (branchData(i,2) == busNumber)
            yShunt = 0.5*branchData(i,6);
            Bi0 = Bi0 + yShunt;
        end
    end

end

%% Group Buses base on Bus Type
function [PQbus,PVbus,SLACKbus] = GroupBuses(busData)
    
    PQbus = busData((busData(:,2) == 0),1);         
    PVbus = busData((busData(:,2) == 2),1);
    SLACKbus = busData((busData(:,2) == 3),1);
    
end
