%%
% FREE CALCIUM CONCENTRATION IN NON-EQUILIBRIUM CONDITIONS
%
% - discretization of the reaction-diffusion equations
% - Ca, CaB and B concentrations temporal and spatial evolution
% - Fluorescence estimation (and its temporal derivative)
% - Comparison of the model for the calcium concentration estimation
%   starting from fluorescence measurements
%
%
%%

clc; clear; close all;
disp('START SIMULATION');

% time evolution:
% - if true: show time-evolution of calcium and free dye spatial concentration
% - if false: hide time-evolution of calcium and free dye spatial concentration
timeEvolution = false;

tic;                                    % starting computational time

%% Physiological constants
Dca = 0.44;                             % calcium diffusion coefficient (microm^2 / ms)

% BAPTA
Db = 0.27;                              % dye diffusion coefficient (microm^2 / ms)
Kplus = 0.5;                            % dye binding rate constant (1 / microM / ms)
Kminus = 0.096;                         % dye unbinding rate constant (1 / ms)
Kd = Kminus/Kplus;                      % dye dissociation constant (microM)

%% Simulation variables initialization
Nvox = 101;                             % total number of voxels
vox = 1:Nvox;
dx = 1.0;                               % voxel size (microm)
dt = 10^-4;                             % integration time (ms)
totalTime = 80;                         % total simulation time (ms)
steps = totalTime / dt;
fluxDuration = 10;                      % total influx time (ms)
fluxEnd = fluxDuration / dt;

J = 50.0;                               % influx current density (microM / ms)
totalJ = J * fluxDuration;              % total entering calcium concentration (microM)
n0 = ceil(Nvox/2);                      % influx position voxel

%% Concentration initialization
Cai = 0.1;                              % initial free calcium concentration (microM)
Btot = 100;                             % total dye concentration (microM)

% concentration initialization for B and CaB -> steady state solution
Bi = Kd / (Kd + Cai) * Btot;            % initial free dye concentration (microM)
CaBi = Btot - Bi;                       % initial bound dye concentration (microM)

% array to store run-time concentrations
Ca = zeros(1,Nvox) + Cai;
B = zeros(1,Nvox)+ Bi;
CaB = zeros(1,Nvox) + CaBi;


%% Discrete reaction-diffusion model
% diffusion probability
pCa = 2 * Dca / dx^2 * dt;              % fraction of calcium transferred to adjacent positions
pB = 2 * Db / dx^2 * dt;                % fraction of dye transferred to adjacent positions

% kernel used to evaluate the laplacian
kernelCa = [pCa/2 -pCa pCa/2];
kernelB = [pB/2 -pB pB/2];

% array to store time evolution concentrations at fixed distance
% for istance, in flux entering position
Cat = zeros(1,steps);                   % temporal evolution of Ca at fixed distance
Cat(1) = Cai;                           % initial condition
CaBt = zeros(1,steps);                  % temporal evolution of CaB at fixed distance
CaBt(1) = CaBi;                         % initial condition
Bt = zeros(1,steps);                    % temporal evolution of B at fixed distance
Bt(1) = Bi;                             % initial condition

% array to store the laplacian of CaB
% used in the fluorescence formula evaluation
laplCaB = zeros(1,steps);

% create empty figures
if timeEvolution == true
    figure('Name','Concentration time evolution','NumberTitle','off','Position',[1 400 600 400]);
    ax1 = subplot(2,1,1);
    grid(ax1,'on');
    ax2 = subplot(2,1,2);
    grid(ax2,'on');
end

% LOOP OVER TIME STEPS
for k = 2:steps

    % Reaction-diffusion for Ca
    laplacianCa = conv2(Ca, kernelCa, 'same');
    
    % reflective boundary
    laplacianCa(1) = -pCa*Ca(1) + pCa*Ca(2);
    laplacianCa(Nvox) = -pCa*Ca(Nvox) + pCa*Ca(Nvox-1);

    R = Kminus .* CaB - Kplus .* Ca .* B;        % reaction term
    
    deltaCa = laplacianCa + R * dt;
    
    if k < fluxEnd
        deltaCa(n0) = deltaCa(n0) + J * dt;
    end
    
    
    % Reaction-diffusion for B
    laplacianB = conv2(B, kernelB, 'same');
    
    % reflective boundary
    laplacianB(1) = -pB*B(1) + pB*B(2);
    laplacianB(Nvox) = -pB*B(Nvox) + pB*B(Nvox-1);

    deltaB = laplacianB + R * dt;
    
    
    % Reaction-diffusion for CaB
    laplacianCaB = conv2(CaB, kernelB, 'same');
    
    % reflective boundary
    laplacianCaB(1) = -pB*CaB(1) + pB*CaB(2);
    laplacianCaB(Nvox) = -pB*CaB(Nvox) + pB*CaB(Nvox-1);
    
    deltaCaB = laplacianCaB -R * dt;
    
    
    % store laplacian of CaB
    laplCaB(k) = laplacianCaB(n0);

    
    % update concentrations
    Ca = Ca + deltaCa;
    B = B + deltaB;
    CaB = CaB + deltaCaB;
    
    
    % Plot
    if timeEvolution == true
        if (mod(k,10^4)==0) && (k*dt < 30)
            plot(ax1, vox, Ca);
            grid(ax1,'on');
            ylim(ax1, [0 100]);
            xlabel(ax1, 'X (\mum)');
            ylabel(ax1, '[Ca] (\muM)')
            if k < fluxEnd
                title(ax1, ['Ca concentration               Time: ',num2str(k*dt),'ms, Flux: ON']);
            else
                title(ax1, ['Ca concentration               Time: ',num2str(k*dt),'ms, Flux: OFF']);
            end

            plot(ax2, vox, CaB);
            grid(ax2,'on');
            ylim(ax2, [20 110]);
            xlabel(ax2, 'X (\mum)');
            ylabel(ax2, '[CaB] (\muM)')
            if k < fluxEnd
                title(ax2, ['CaB concentration              Time: ',num2str(k*dt),'ms, Flux: ON']);
            else
                title(ax2, ['CaB concentration              Time: ',num2str(k*dt),'ms, Flux: OFF']);
            end
            pause(0.4);
        end
    end
    
    % evaluate the constancy of total dye concentration in each position
    %if mod(k,500)==0
    %    temp = B + CaB
    %end
    
    % update temporal evolution in flux enter position
    Cat(k) = Ca(n0);
    CaBt(k) = CaB(n0);
    Bt(k) = B(n0);
end
% END LOOP OVER TIME STEPS


%% Concentration at fixed distance
figure('Name','Concentration at fixed distance','NumberTitle','off','Position',[800 400 700 400]);
ax3 = subplot(2,1,1);
ax4 = subplot(2,1,2);

% Free calcium concentration plot
plot(ax3, (1:steps)*dt, Cat, '-b');
xlabel(ax3, 'Time (ms)');
ylabel(ax3, '[Ca] (\muM)');
title(ax3, 'Ca concentration at flux entering position');
xline(ax3,(fluxEnd-0.5)*dt,'--');
ylim(ax3, [min(Cat)-max(Cat)/10 max(Cat)+max(Cat)/10]);

% Calcium bound concentration plot
plot(ax4, (1:steps)*dt, CaBt, '-b');
xlabel(ax4, 'Time (ms)');
ylabel(ax4, '[CaB] (\muM)');
title(ax4, 'CaB concentration at flux entering position');
xline(ax4,(fluxEnd-0.5)*dt,'--');
ylim(ax4, [min(CaBt)-max(CaBt)/10 max(CaBt)+max(CaBt)/10]);


%% Fluorescence simulation
Sf = 1;
Sb = 5;

% Fluorescence
Fmax = Sb * Btot;
Fmin = Sf * Btot;
F = Sf * Bt + Sb * CaBt;

% Fluorescence temporal derivative
kernelDer = [1 -1] / dt;
dF = conv2(F, kernelDer, 'same');   % discretization of temporal derivative
dFi = 0.0;                          % since initial configuration is in equilibrium
dF = [dFi dF];
dF(end) = [];

figure('Name','Fluorescence','NumberTitle','off','Position',[1 400 700 400]);
ax5 = subplot(2,1,1);
ax6 = subplot(2,1,2);

% Fluorescence signal plot
plot(ax5,(1:steps)*dt, F, '-b', 1:totalTime, Fmin+zeros(1,totalTime), '-.k', 1:totalTime, Fmax+zeros(1,totalTime), '--k');
xlabel(ax5, 'Time (ms)');
ylabel(ax5, 'F (A.U.)');
ylim(ax5, [Fmin-80 Fmax+80]);
title(ax5, 'Fluorescence');
legend(ax5, 'F', 'Fmin', 'Fmax');

% Fluorescence temporal derivative plot
plot(ax6,(1:steps)*dt, dF, '-b');
xlabel(ax6, 'Time (ms)');
ylabel(ax6, 'dF/dt (A.U./ms)');
title(ax6, 'Fluorescence derivative');


%% Models evaluation

% Calcium concentration from estimated formulae
CaEq = Kd .* (F-Fmin) ./ (Fmax - F);
CaNoLapl = (dF + (F-Fmin).*Kminus) ./ Kplus ./ (Fmax-F);
A = (Fmax - Fmin)/Btot;
CaLapl = (dF - laplCaB .* A ./ dt + (F-Fmin).*Kminus) ./ Kplus ./ (Fmax-F);

% Models comparing plot 
figure('Name','Model evaluation','NumberTitle','off','Position', [400 1 700 400]);
plot((1:steps)*dt, CaEq, '-b', (1:steps)*dt, CaNoLapl, '-m', (1:steps)*dt, CaLapl, '-r', 'LineWidth',1);
hold on;
plot((1:5000:fluxEnd+15000)*dt, Cat(1:5000:fluxEnd+15000), 'Linestyle', 'none', 'Marker', '.', 'MarkerSize',10,'MarkerEdgeColor','k');
plot((fluxEnd+15000:20000:steps)*dt, Cat(fluxEnd+15000:20000:steps), 'Linestyle', 'none', 'Marker', '.', 'MarkerSize',10,'MarkerEdgeColor','k');
ylim([min(Cat)-max(Cat)/10 max(Cat)+max(Cat)/10]);
xlabel('Time (ms)');
ylabel('[Ca] (\muM)');
title('Model evaluation');
legend('Eq formula','Non-eq formula wo lapl','Non-eq formula w lapl','[Ca](t)');


%% Evaluate computation time
computationTime = toc;
disp('END SIMULATION');
disp(['Total computation time: ', num2str(computationTime),' s']);