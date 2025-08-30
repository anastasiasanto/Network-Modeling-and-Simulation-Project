%MODELLO AGGREGATION ONLY
% rateConstants are the rate constants for each reaction.
% stoichiometry is a 2D array where each row corresponds to a reaction, and each column
% corresponds to a species. The values represent the change in the number of molecules of each species due to the reaction.
% initialState is the initial number of molecules for each species.
% timeEnd is the simulation time.
% The function returns the simulation results: times is the time points, states is the number of molecules
% for each species at each time point, reactionsOccurred is the index of the reaction that occurred at each step,
% and reactionTimes is the time at which each reaction occurred. The states matrix is transposed to match the pseudocode.

Stoichiometrynofrag = [
    -2, 1, 0, 0, 0, 0, 0, 0;
    -1, -1, 1, 0, 0, 0, 0, 0;
    -1, 0, -1, 1, 0, 0, 0, 0;
    0, -2, 0, 1, 0, 0, 0, 0;
    -1, 0, 0, -1, 1, 0, 0, 0;
    -1, -2, 0, 0, 1, 0, 0, 0;
    0, -1, -1, 0, 1, 0, 0, 0;
    -1, 0, 0, 0, -1, 1, 0, 0;
    0, -1, 0, -1, 0, 1, 0, 0;
    0, 0, -2, 0, 0, 1, 0, 0;
    -2, -2, 0, 0, 0, 1, 0, 0;
    -1, 0, 0, 0, 0, -1, 1, 0;
    0, -1, 0, 0, -1, 0, 1, 0;
    -1, 0, -2, 0, 0, 0, 1, 0;
    0, 0, -1, -1, 0, 0, 1, 0;
    -1, -3, 0, 0, 0, 0, 1, 0;
    -1, 0, 0, 0, 0, 0, -1, 1;
    0, -1, 0, 0, 0, -1, 0, 1;
    0, 0, -1, 0, -1, 0, 0, 1;
    0, 0, 0, -2, 0, 0, 0, 1;
    -1, -1, 0, 0, -1, 0, 0, 1;
    -1, 0, -1, -1, 0, 0, 0, 1;
    0, -2, 0, -1, 0, 0, 0, 1;
    0, -1, -2, 0, 0, 0, 0, 1;
    -2, 0, -2, 0, 0, 0, 0, 1;
    0, -4, 0, 0, 0, 0, 0, 1;
    -2, -1, 0, -1, 0, 0, 0, 1
];

% The table S1 only contains the kinetics of amyloid-beta aggregation and
% does not consider the fragmentation.
% M1 + M1 = M2 : k0M1M1
% M2 + M2 + M1 + M1 = M6 : k10M2M2M1M1
% M1 + M2 = M3 :k1M1M2
% M1 + M6 = M7 : k11M1M6
% M1 + M3 = M4 : k2M1M3
% M2 + M5 = M7 : k12M2M5
% M2 + M2 = M4 : k3M2M2
% M3 + M3 + M1 = M7 : k13M3M3M1
% M1 + M4 = M5 : k4M1M4
% M3 + M4 = M7 : k14M3M4
% M1 + M2 + M2 = M5 : k5M1M2M2
% M1 + 3M2 = M7 : k15M1M2M2M2
% M3 + M3 = M5 : k6M2M3
% M1 + M7 = M8 : k16M1M7
% M1 + M5 = M6 : k7M1M5
% M2 + M6 = M8 : k17M2M6
% M2 + M4 = M6 : k8M2M6
% M3 + M5 = M8 : k18M3M5
% M3 + M3 = M6 : k9M3M3
% M4 + M4 = M8 : k19M4M4
% M1 + M2 + M5 = M8 : k20M1M2M5
% M1 + M3 + M4 = M8 : k21M1M3M4
% M2 + M2 + M4 = M8 : k22M2M2M4
% M2 + M3 + M3 = M8 : k23M2M3M3
% M1 + M1 + M3 + M3 = M8 : k24M1M1M3M3
% 4M2 = M8 : k25M2M2M2M2
% M1 + M1 + M2 + M4 = M8 : k26M1M1M2M

%in figura 3 abbiamo rappresentato due condizioni in cui : m1=10000, m1=
%5000; the initial condition for the remaining population m2,...,m8 -->
%100,60,50,20,10,5,5. k0=0.00001, k1=1/2k0 e così via...

rateconstants = [1e-05; 5e-06; 2.5e-06; 1.25e-06; 6.25e-07;3.125e-07; 1.5625e-07; 7.8125e-08; 3.90625e-08; 1.953125e-08;
    9.765625e-09; 4.8828125e-09; 2.44140625e-09; 1.220703125e-09; 6.103515625e-10; 3.0517578125e-10; 1.52587890625e-10; 
    7.62939453125e-11; 3.814697265625e-11; 1.9073486328125e-11; 9.5367431640625e-12; 4.76837158203125e-12; 2.384185791015625e-12;
    1.1920928955078126e-12; 5.960464477539063e-13; 2.9802322387695315e-13; 1.4901161193847657e-13];

InitialConditions31=[10000,100,60,50,20,10,5,5];%Condition1
InitialConditions32=[5000,100,60,50,20,10,5,5];%Condition2
% Inizializzazione della prima k
k0 = 0.00001;
% k1 = 5e-06;
% k2 = 2.5e-06;
% k3 = 1.25e-06;
% k4 = 6.25e-07;
% k5 = 3.125e-07;
% k6 = 1.5625e-07;
% k7 = 7.8125e-08;
% k8 = 3.90625e-08;
% k9 = 1.953125e-08;
% k10 = 9.765625e-09;
% k11 = 4.8828125e-09;
% k12 = 2.44140625e-09;
% k13 = 1.220703125e-09;
% k14 = 6.103515625e-10;
% k15 = 3.0517578125e-10;
% k16 = 1.52587890625e-10;
% k17 = 7.62939453125e-11;
% k18 = 3.814697265625e-11;
% k19 = 1.9073486328125e-11;
% k20 = 9.5367431640625e-12;
% k21 = 4.76837158203125e-12;
% k22 = 2.384185791015625e-12;
% k23 = 1.1920928955078126e-12;
% k24 = 5.960464477539063e-13;
% k25 = 2.9802322387695315e-13;
% k26 = 1.4901161193847657e-13; 

%in figura 4 cambia solo il fatto che m2=...=m8=1, quindi il resto è uguale
%ma cambiano le condizioni iniziali :
InitialConditions41=[10000,1,1,1,1,1,1,1];%Condition1
InitialConditions42=[5000,1,1,1,1,1,1,1];%Condition2
% Esegui la simulazione
timeEnd=1000;
[timeGrid, avgStates, avgStd, ~, ~] = MonteCarloSimulation(rateconstants, Stoichiometrynofrag,InitialConditions41 , timeEnd,100);
[timeGrid2, avgStates2, avgStd2, ~, ~] = MonteCarloSimulation(rateconstants, Stoichiometrynofrag,InitialConditions42 , timeEnd,100);
% Plotta i risultati medi e ombreggia la deviazione standard
% Creare la prima figura con 4 subplot
figure(1);

for i = 1:4
    subplot(2, 2, i);
    
    % Plotta la dinamica media
    plot(timeGrid, squeeze(avgStates(1, i, :)), 'b', 'LineWidth', 2);
    hold on;
    
    % Calcola i limiti superiori e inferiori dell'area ombreggiata
    upperBound = squeeze(avgStates(1, i, :) + avgStd(1, i, :));
    lowerBound = squeeze(avgStates(1, i, :) - avgStd(1, i, :));

    % Regola la lunghezza dei vettori per evitare errori
    timeGridAdjusted = linspace(0, timeEnd, length(upperBound));

    % Assicurati che tutti i vettori abbiano la stessa lunghezza
    upperBound = upperBound(:)';
    lowerBound = lowerBound(:)';

    % Ombreggia l'area tra la deviazione standard superiore e inferiore
    fill([timeGridAdjusted, fliplr(timeGridAdjusted)], [upperBound, fliplr(lowerBound)], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    
    title(['Species ' num2str(i)]);
    xlabel('Time');
    ylabel('Population');
    legend('Mean', 'Standard Deviation', 'Location', 'Best');
end

% Creare la seconda figura con i rimanenti 4 subplot
figure(2);

for i = 5:size(avgStates, 2)
    subplot(2, 2, i-4);
    
    % Plotta la dinamica media
    plot(timeGrid, squeeze(avgStates(1, i, :)), 'b', 'LineWidth', 2);
    hold on;
    
    % Calcola i limiti superiori e inferiori dell'area ombreggiata
    upperBound = squeeze(avgStates(1, i, :) + avgStd(1, i, :));
    lowerBound = squeeze(avgStates(1, i, :) - avgStd(1, i, :));

    % Regola la lunghezza dei vettori per evitare errori
    timeGridAdjusted = linspace(0, timeEnd, length(upperBound));

    % Assicurati che tutti i vettori abbiano la stessa lunghezza
    upperBound = upperBound(:)';
    lowerBound = lowerBound(:)';

    % Ombreggia l'area tra la deviazione standard superiore e inferiore
    fill([timeGridAdjusted, fliplr(timeGridAdjusted)], [upperBound, fliplr(lowerBound)], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    
    title(['Species ' num2str(i)]);
    xlabel('Time');
    ylabel('Population');
    legend('Mean', 'Standard Deviation', 'Location', 'Best');
end

% Aggiungi i risultati della seconda simulazione alle figure esistenti
figure(1); % Seleziona la prima figura
for i = 1:4
    subplot(2, 2, i);
    hold on;
    % Plotta la dinamica media della seconda simulazione
    plot(timeGrid2, squeeze(avgStates2(1, i, :)), 'g', 'LineWidth', 2); % Usa il colore verde per distinguere la seconda simulazione
    legend('Mean (Condition 1)', 'Std Dev (Condition 1)', 'Mean (Condition 2)', 'Location', 'Best');
    upperBound = squeeze(avgStates2(1, i, :) + avgStd2(1, i, :));
    lowerBound = squeeze(avgStates2(1, i, :) - avgStd2(1, i, :));

    % Regola la lunghezza dei vettori per evitare errori
    timeGridAdjusted = linspace(0, timeEnd, length(upperBound));

    % Assicurati che tutti i vettori abbiano la stessa lunghezza
    upperBound = upperBound(:)';
    lowerBound = lowerBound(:)';

    % Ombreggia l'area tra la deviazione standard superiore e inferiore
    fill([timeGridAdjusted, fliplr(timeGridAdjusted)], [upperBound, fliplr(lowerBound)], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
   
end

figure(2); % Seleziona la seconda figura
for i = 5:size(avgStates2, 2)
    subplot(2, 2, i-4);
    hold on;
    % Plotta la dinamica media della seconda simulazione
    plot(timeGrid2, squeeze(avgStates2(1, i, :)), 'g', 'LineWidth', 2); % Usa il colore verde per distinguere la seconda simulazione
    legend('Mean (Condition 1)', 'Std Dev (Condition 1)', 'Mean (Condition 2)', 'Location', 'Best');
    upperBound = squeeze(avgStates2(1, i, :) + avgStd2(1, i, :));
    lowerBound = squeeze(avgStates2(1, i, :) - avgStd2(1, i, :));

    % Regola la lunghezza dei vettori per evitare errori
    timeGridAdjusted = linspace(0, timeEnd, length(upperBound));

    % Assicurati che tutti i vettori abbiano la stessa lunghezza
    upperBound = upperBound(:)';
    lowerBound = lowerBound(:)';

    % Ombreggia l'area tra la deviazione standard superiore e inferiore
    fill([timeGridAdjusted, fliplr(timeGridAdjusted)], [upperBound, fliplr(lowerBound)], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
   

end

% Controlla le dimensioni di  avgStates e  timeGrid
[numSimulations, numSpecies, numTimePoints] = size(avgStates);
if numel(timeGrid) ~= numTimePoints
    error('Il numero di time points non fitta la misura di timeGrid.');
end

% Inizializza i vettori per il momento zero e il primo momento
zerothMoment = zeros(1, size(avgStates2, 3));
firstMoment = zeros(1, size(avgStates2, 3));

% Calcola il momento zero e il primo momento per ogni istante di tempo
for t = 1:size(avgStates2, 3)
    zerothMoment(t) = sum(avgStates2(1, 2:end, t));
    for s = 2:size(avgStates2, 2) % Itera da 2 a size(avgStates, 2) per escludere la prima specie
        firstMoment(t) = firstMoment(t) + sum(avgStates2(1, s, t) * s);
    end
end

% Crea una nuova figura
figure;

% Plotta il momento zero
subplot(2, 1, 1);
plot(1:size(avgStates2, 3), zerothMoment, 'b', 'LineWidth', 2);
title('Zeroth Moment');
xlabel('Time');
ylabel('Total Population (excluding species 1)');

% Plotta il primo momento
subplot(2, 1, 2);
plot(1:size(avgStates2, 3), firstMoment, 'r', 'LineWidth', 2);
title('First Moment');
xlabel('Time');
ylabel('First Moment (Mass Concentration) (excluding species 1)');
