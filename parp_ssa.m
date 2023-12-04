% Simulate a two-state model of gene expression
import Gillespie.*

%% Species network:
%   1. Free DNA
%   2. Free Prot
%   3. DNA Prot Complex:
%   4*(Round-1)+3. DNA-PARP Round
%   4*(Round-1)+4. DNA-PARP-PARPi Round
%   4*(Round-1)+5. PARP Round
%   4*(Round-1)+6. PARP-PARPi Round
%   End. Final PARP
%% Rate constants
iterations = 100;
end_time = 7200;

%Things that can be changed
drugConcentration = 1e-6; % ~1 uM saturates target ::: CHECK with initial values

%fitted through modeling to match 80% recovery @ 3-10 minutes
rounds = 50;%50
p.PARPo = 5; % ~1uM %5
p.Proto = 1; % protein concentration %1 set to 0.5 in iteration loop
p.gamma = 2; %lower PARP1-DNA binding upon parylation %2
p.prot_thresh = 2; %2
p.Histones = 10*p.PARPo; % other proteins that can be parylated %10*p.PARPo

%Constants
p.DNAo = 1;
p.Drug = drugConcentration*1e6*p.PARPo; % relate to p.PARPo
p.NAD = 100*p.PARPo; %100 uM
p.rounds = rounds;
 
%Rates - ka is 10-7 * actual rate
%drug constants - independent of PAR
p.kd_trap = 8e-5; % kd trap   vel - 5e-3; olap - 3e-4; tal - 8e-5;
p.ka_trap = 3.6e-2; % ka trap vel - 1.8e-1; olap - 2e-2; tal - 3.6e-2;

%NAD constants
p.ka_PAR = 5e-2;%5e-2

%Other protein binding
p.ka_Prot = 1e-2; % protein binding DNA 1e-2
p.kd_Prot = 4e-4; %PARP releasing from DNA one order of magnitude stronger 4e-4
 
%Initial PARP DNA dissociation
p.kd_rel0 = 4e-3; %PARP releasing from DNA
p.ka_rel0 = 1e4; % PARP binding DNA
 
p.kd_relD0 = 2.52e-3; %PARP-PARPi releasing from DNA % zandarashvili vel - 4.08e-3; olap - 2.31e-3; tal - 2.52e-3
p.ka_relD0 = 1e-2; % PARP-PARPi binding DNA
 
%Exchange rate
p.exchange = 0.2/p.PARPo; %20% takes ~2 sec Shao 9694?9709 NAR, 2020, Vol. 48, No. 17

 
%% Initial state
tspan = [0, end_time];
x0    = zeros(1,4+rounds*4);
x0(1) = p.DNAo;
x0(2) = p.Proto;
x0(6) = p.PARPo; % set based on amount of drug in system 6 for no drug in system 7 for drug
%% Specify reaction network
pfun = @propensities_2state;

p.react_number = 11;
stoich_matrix = spalloc(4+rounds*p.react_number, 4+rounds*4,3*(4+rounds*p.react_number));
stoich_matrix(1,[1,2,3]) = [-1,-1,1]; % DNA protein biding 1 + 2 --> 3
stoich_matrix(2,[1,2,3]) = [1,1,-1]; % dissociation of DNA protein 3 --> 1+2
for i = 1:rounds
    stoich_matrix(p.react_number*(i-1)+3,[1,4*(i-1)+4,4*(i-1)+6]) = [-1,1,-1];% 1 + 6 --> 4
    stoich_matrix(p.react_number*(i-1)+4,[1,4*(i-1)+4,4*(i-1)+6]) = [1,-1,1];% 4 --> 1 + 6
    stoich_matrix(p.react_number*(i-1)+5,[4*(i-1)+4,4*(i-1)+5]) = [-1,1];% 4 --> 5
    stoich_matrix(p.react_number*(i-1)+6,[4*(i-1)+4,4*(i-1)+5]) = [1,-1];% 5 --> 4
    stoich_matrix(p.react_number*(i-1)+7,[1,4*(i-1)+5,4*(i-1)+7]) = [1,-1,1];% 5 --> 1 + 7
    stoich_matrix(p.react_number*(i-1)+8,[1,4*(i-1)+5,4*(i-1)+7]) = [-1,1,-1];% 1 + 7 --> 5
    stoich_matrix(p.react_number*(i-1)+9,:) = 0; % PAR
    stoich_matrix(p.react_number*(i-1)+10, [4*(i-1)+6,4*(i-1)+7]) = [-1,1]; % 6 --> 7
    stoich_matrix(p.react_number*(i-1)+11,[4*(i-1)+6,4*(i-1)+7]) = [1,-1]; % 7 --> 6
    stoich_matrix(p.react_number*(i-1)+12,[6,(4*(i-1)+6)]) = [1,-1]; % 6 parp exchange
    stoich_matrix(p.react_number*(i-1)+13,[7,(4*(i-1)+7)]) = [1,-1]; % 7 parp exchange
end
stoich_matrix(end-1,[6,end]) = [1,-1]; % 6 parp exchange with saturated
stoich_matrix(end,[1,end-4,end]) = [1,-1,1]; % pentultimate --> 1 and Saturated 

%% Run simulation
final_x = zeros(tspan(end),4+rounds*4);
final_par_count = zeros(tspan(end),1);
final_system_par = zeros(tspan(end),1);
%final_prot_count = zeros(tspan(end),1);

stop_time = tspan(end)*ones(iterations,1);
exchange_counter = zeros(iterations,1);

prot_count_iteration = zeros(iterations,tspan(end));
par_count_iteration = zeros(iterations,tspan(end));

for i = 1:iterations
    if i <= iterations/2 %set average ptot0 to 0.5 - don't want fractions
        x0(2) = 0;
    else
        x0(2) = p.Proto;
    end
    prelim_x = zeros(tspan(end),4+rounds*4); % added to store data to add to final
    prelim_par_count = zeros(tspan(end),1);
    prelim_system_par = zeros(tspan(end),1);
    %prelim_prot_count = zeros(tspan(end),1);
    
    [t,x, par_counter, par_in_system, exchange_counter(i)] = directMethod(stoich_matrix, pfun, tspan, x0, p);
    
    new_end = ceil(min([tspan(end),stop_time(i)]));
    
    time_bins = discretize(t,new_end);
    unique_time = unique(time_bins);
    
    %minimum of tspan/stop time?
    for j = 1:tspan(end)
        if j == 1
            prelim_x(j,:) = x0; 
            prelim_par_count(1) = par_counter(1);
            prelim_system_par(1) = par_in_system(1);
            %prelim_prot_count(1) = prot_counter(1);
            
        elseif ismember(j,unique_time) == 1
            rowValue = find(time_bins == j,1, 'last'); %pick out last value in tstep
            prelim_x(j,:) = x(rowValue,:);
            prelim_par_count(j) = par_counter(rowValue);
            prelim_system_par(j) = par_in_system(rowValue);
            %prelim_prot_count(j) = prot_counter(rowValue);
            
        else
            prelim_x(j,:) = prelim_x((j-1),:);
            prelim_par_count(j) = prelim_par_count(j-1);
            prelim_system_par(j) = prelim_system_par(j-1);
            %prelim_prot_count(j) = prelim_prot_count(j-1);
        end
    end

    final_x = final_x + prelim_x;
    final_par_count = final_par_count + prelim_par_count;
    final_system_par = final_system_par + prelim_system_par;
    %final_prot_count = final_prot_count + prelim_prot_count;
    
    prot_count_iteration(i,:) = prelim_x(:,2);
    par_count_iteration(i,:) = prelim_system_par(:,1);
end

%% Things to care about
final_x = final_x/iterations;
final_par_count = final_par_count/iterations;
final_system_par = final_system_par/iterations;


final_trapped = sum(final_x(:,5:4:end-3),2);
final_parp = sum(final_x(:,4:4:end-4),2);
final_parpFree = sum(final_x(:,6:4:end-2),2);
final_parpDrug = sum(final_x(:,7:4:end-1),2);
final_totParp = final_trapped+final_parp+final_parpFree+final_parpDrug+final_x(:,end);
final_dnaTotal = final_x(:,1)+final_trapped+final_parp;
final_parpOnDNA = final_parp+final_trapped;
final_prot_count = final_x(:,2)+final_x(:,3);%final_prot_count/iterations;

round_onDNA = zeros(end_time,rounds);
round_free = zeros(end_time,rounds);

% plotting heatmap
for z=1:rounds
    round_onDNA(:,z) = final_x(:,(4+4*(z-1)))+final_x(:,(5+4*(z-1)));
    round_free(:,z) = final_x(:,(6+4*(z-1)))+final_x(:,(7+4*(z-1)));
end

round_parpTot = round_onDNA + round_free;
round_parpPAR = round_parpTot(:,2:50);
    

%% Plotting
figure(1);
hold on
plot(1:tspan(end),final_par_count,'r'); plot(1:tspan(end),final_system_par,'k'); 
legend('parCount','systemPar');
hold off
 
figure(2);
%colormap('hot')
imagesc(round_parpTot)
xlabel('species')
ylabel('time (s)')
colorbar
 
figure(3);
imagesc(round_parpPAR)
xlabel('species')
ylabel('time (s)')
colorbar



%figure(3);
%hold on
%plot(1:tspan(end),final_x(:,end))
%legend('saturated PARP');
%hold off
 
figure(4);
hold on
plot(1:tspan(end),final_parpOnDNA,'k');
legend('onDNA');
hold off
 
figure(5);
hold on
plot(1:tspan(end),final_prot_count,'g'); 
legend('protein');
hold off

%% Helper functions
function a = propensities_2state(x, p)
% Return reaction propensities given current state x

p.kd_rel_logscale = logspace(log10(p.kd_rel0),log10(p.kd_rel0)+p.gamma,p.rounds);
p.ka_rel_logscale = logspace(log10(p.ka_rel0),log10(p.ka_rel0)-p.gamma,p.rounds);

p.kd_relD_logscale = logspace(log10(p.kd_relD0),log10(p.kd_relD0)+p.gamma,p.rounds);
p.ka_relD_logscale = logspace(log10(p.ka_relD0),log10(p.ka_relD0)-p.gamma,p.rounds);

a = zeros(4+p.rounds*p.react_number,1);
a(1)= p.ka_Prot*x(1)*x(2);
a(2)=p.kd_Prot*x(3);

for i = 1:p.rounds
    a(p.react_number*(i-1)+3) = p.ka_rel_logscale(i)*x(1)*x(4*(i-1)+6);
    a(p.react_number*(i-1)+4) = p.kd_rel_logscale(i)*x(4*(i-1)+4);
    a(p.react_number*(i-1)+5) = p.ka_trap*x(4*(i-1)+4)*p.Drug;
    a(p.react_number*(i-1)+6) = p.kd_trap*x(4*(i-1)+5);
    a(p.react_number*(i-1)+7) = p.kd_relD_logscale(i)*x(4*(i-1)+5);
    a(p.react_number*(i-1)+8) = p.ka_relD_logscale(i)*x(4*(i-1)+7)*x(1);
    a(p.react_number*(i-1)+9) = p.ka_PAR*p.NAD*x(4*(i-1)+4);
    a(p.react_number*(i-1)+10) = p.ka_trap*x(4*(i-1)+6)*p.Drug;
    a(p.react_number*(i-1)+11) = p.kd_trap*x(4*(i-1)+7);
    a(p.react_number*(i-1)+12) = p.exchange*x(4*(i-1)+6);
    a(p.react_number*(i-1)+13) = p.exchange*x(4*(i-1)+7);
end

a(end-1)= p.exchange*x(end);
a(end) = p.ka_PAR*p.NAD*x(end-4);

a(12:13) = 0;

end