%MODELLO AGGREGATION-FRAGMENMTATION:
% rateConstants are the rate constants for each reaction.
% stoichiometry is a 2D array where each row corresponds to a reaction, and each column
% corresponds to a species. The values represent the change in the number of molecules of each species due to the reaction.
% initialState is the initial number of molecules for each species.
% timeEnd is the simulation time.
% The function returns the simulation results: times is the time points, states is the number of molecules
% for each species at each time point, reactionsOccurred is the index of the reaction that occurred at each step,
% and reactionTimes is the time at which each reaction occurred. The states matrix is transposed to match the pseudocode.

Stoichiometry2=[
    -2.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00;
    2.00, -1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00;
    -1.00, -1.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00;
    1.00, 1.00, -1.00, 0.00, 0.00, 0.00, 0.00, 0.00;
    -1.00, 0.00, -1.00, 1.00, 0.00, 0.00, 0.00, 0.00;
    1.00, 0.00, 1.00, -1.00, 0.00, 0.00, 0.00, 0.00;
    0.00, -2.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00;
    0.00, 2.00, 0.00, -1.00, 0.00, 0.00, 0.00, 0.00;
    -1.00, 0.00, 0.00, -1.00, 1.00, 0.00, 0.00, 0.00;
    1.00, 0.00, 0.00, 1.00, -1.00, 0.00, 0.00, 0.00;
    -1.00, -2.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00;
    1.00, 2.00, 0.00, 0.00, -1.00, 0.00, 0.00, 0.00;
    0.00, -1.00, -1.00, 0.00, 1.00, 0.00, 0.00, 0.00;
    0.00, 1.00, 1.00, 0.00, -1.00, 0.00, 0.00, 0.00;
    -1.00, 0.00, 0.00, 0.00, -1.00, 1.00, 0.00, 0.00;
    1.00, 0.00, 0.00, 0.00, 1.00, -1.00, 0.00, 0.00;
    0.00, -1.00, 0.00, -1.00, 0.00, 1.00, 0.00, 0.00;
    0.00, 1.00, 0.00, 1.00, 0.00, -1.00, 0.00, 0.00;
    0.00, 0.00, -2.00, 0.00, 0.00, 1.00, 0.00, 0.00;
    0.00, 0.00, 2.00, 0.00, 0.00, -1.00, 0.00, 0.00;
    -2.00, -2.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00;
    2.00, 2.00, 0.00, 0.00, 0.00, -1.00, 0.00, 0.00;
    -1.00, 0.00, 0.00, 0.00, 0.00, -1.00, 1.00, 0.00;
    1.00, 0.00, 0.00, 0.00, 0.00, 1.00, -1.00, 0.00;
    0.00, -1.00, 0.00, 0.00, -1.00, 0.00, 1.00, 0.00;
    0.00, 1.00, 0.00, 0.00, 1.00, 0.00, -1.00, 0.00;
    -1.00, 0.00, -2.00, 0.00, 0.00, 0.00, 1.00, 0.00;
    1.00, 0.00, 2.00, 0.00, 0.00, 0.00, -1.00, 0.00;
    0.00, 0.00, -1.00, -1.00, 0.00, 0.00, 1.00, 0.00;
    0.00, 0.00, 1.00, 1.00, 0.00, 0.00, -1.00, 0.00;
    -1.00, -3.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00;
    1.00, 3.00, 0.00, 0.00, 0.00, 0.00, -1.00, 0.00;
    -1.00, 0.00, 0.00, 0.00, 0.00, 0.00, -1.00, 1.00;
    1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, -1.00;
    0.00, -1.00, 0.00, 0.00, 0.00, -1.00, 0.00, 1.00;
    0.00, 1.00, 0.00, 0.00, 0.00, 1.00, 0.00, -1.00;
    0.00, 0.00, -1.00, 0.00, -1.00, 0.00, 0.00, 1.00;
    0.00, 0.00, 1.00, 0.00, 1.00, 0.00, 0.00, -1.00;
    0.00, 0.00, 0.00, -2.00, 0.00, 0.00, 0.00, 1.00;
    0.00, 0.00, 0.00, 2.00, 0.00, 0.00, 0.00, -1.00;
    -1.00, -1.00, 0.00, 0.00, -1.00, 0.00, 0.00, 1.00;
    1.00, 1.00, 0.00, 0.00, 1.00, 0.00, 0.00, -1.00;
    -1.00, 0.00, -1.00, -1.00, 0.00, 0.00, 0.00, 1.00;
    1.00, 0.00, 1.00, 1.00, 0.00, 0.00, 0.00, -1.00;
    0.00, -2.00, 0.00, -1.00, 0.00, 0.00, 0.00, 1.00;
    0.00, 2.00, 0.00, 1.00, 0.00, 0.00, 0.00, -1.00;
    0.00, -1.00, -2.00, 0.00, 0.00, 0.00, 0.00, 1.00;
    0.00, 1.00, 2.00, 0.00, 0.00, 0.00, 0.00, -1.00;
    -2.00, 0.00, -2.00, 0.00, 0.00, 0.00, 0.00, 1.00;
    2.00, 0.00, 2.00, 0.00, 0.00, 0.00, 0.00, -1.00;
    0.00, -4.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00;
    0.00, 4.00, 0.00, 0.00, 0.00, 0.00, 0.00, -1.00;
    -2.00, -1.00, 0.00, -1.00, 0.00, 0.00, 0.00, 1.00;
    2.00, 1.00, 0.00, 1.00, 0.00, 0.00, 0.00, -1.00
];

% The table S2 only contains the kinetics of amyloid-beta aggregation and
% fragmentation.
% M1 + M1 = M2 : k0M1M1
% M2 = 2M1 :k1M2
% M1 + M2 = M3 :  k2M1M2
% M3 = M1 + M2 : k3M3
% M1 + M3 = M4 : k4M1M3
% M4 = M1 + M3 : k5M4
% M2 + M2 = M4 : k6M2M2
% M4 = 2M2 : k7M4
% M1 + M4 = M5 : k8M1M4
% M5 = M1 + M4 : k9M5
% M1 + M2 + M2 = M5 : k10M1M2M2
% M5 = M1 + 2M2 : k11M5
% M2 + M3 = M5 : k12M2M3
% M5 = M2 + M3 : k13M5
% M1 + M5 = M6 : k14M1M5
% M6 = M1 + M5 : k15M6
% M2 + M4 = M6 : k16M2M4
% M6 = M2 + M4 : k17M6
% M3 + M3 = M6 : k18M3M3
% M6 = 2M3 : k19M6
% M2 + M2 + M1 + M1 = M6 : k20M2M2M1M1
% M6 = 2M2 + 2M1 : k21M6
% M1 + M6 = M7 : k22M1M6
% M7 = M1 + M6 : k23M7
% M2 + M5 = M7 : k24M2M5
% M7 = M2 + M5 : k25M7
% M3 + M3 + M1 = M7 : k26M3M3M1
% M7 = 2M3 + M1 : k27M7
% M3 + M4 = M7 : k28M3M4
% M7 = M3 + M4 : k29M7
% M1 + 3M2 = M7 : k30M1M2M2M2
% M7 = M1 + 3M2 : k31M7
% M1 + M7 = M8 : k32M1M7
% M8 = M1 + M7 : k33M8
% M2 + M6 = M8 : k34M2M6
% M8 = M2 + M6 : k35M8
% M3 + M5 = M8 : k36M3M5
% M8 = M3 + M5 : k37M8
% M4 + M4 = M8 : k38M4M4
% M8 = 2M4 : k39M8
% M1 + M2 + M5 = M8 : k40M1M2M3
% M8 = M1 + M2 + M5 : k41M8
% M1 + M3 + M4 = M8 : k42M1M3M4
% M8 = M1 + M3 + M4 : k43M8
% M2 + M2 + M4 = M8 : k44M2M2M4
% M8 = 2M2 + M4 : k45M8
% M2 + M3 + M3 = M8 : k46M2M3M3
% M8 = M2 + 2M3 : k47M8
% M1 + M1 + M3 + M3 = M8 : k48M1M1M3M3
% M8 = 2M1 + 2M3 : k49M8
% 4M2 = M8 : k50M2M2M2M2
% M8 = 4M2 : k51M8
% M1 + M1 + M2 + M4 = M8: k52M1M1M2M4
% M8 = 2M1 + M2 + M4 : k53M8

%in figura 5 abbiamo rappresentato due condizioni in cui : m1=10000, m1=
%5000; the initial condition for the remaining population m2,...,m8 -->
%100,100,100,100,100,100,100. k0=0.00001, k1=k0/2 e così via...
%RATE CONSTANTS
% k 0 = 0.00001
% k 1 = 5e-06
% k 2 = 5e-06
% k 3 = 2.5e-06
% k 4 = 2.5e-06
% k 5 = 1.25e-06
% k 6 = 1.25e-06
% k 7 = 6.25e-07
% k 8 = 6.25e-07
% k 9 = 3.125e-07
% k 10 = 3.125e-07
% k 11 = 1.5625e-07
% k 12 = 1.5625e-07
% k 13 = 7.8125e-08
% k 14 = 7.8125e-08
% k 15 = 3.90625e-08
% k 16 = 3.90625e-08
% k 17 = 1.953125e-08
% k 18 = 1.953125e-08
% k 19 = 9.765625e-09
% k 20 = 9.765625e-09
% k 21 = 4.8828125e-09
% k 22 = 4.8828125e-09
% k 23 = 2.44140625e-09
% k 24 = 2.44140625e-09
% k 25 = 1.220703125e-09
% k 26 = 1.220703125e-09
% k 27 = 6.103515625e-10
% k 28 = 6.103515625e-10
% k 29 = 3.0517578125e-10
% k 30 = 3.0517578125e-10
% k 31 = 1.52587890625e-10
% k 32 = 1.52587890625e-10
% k 33 = 7.62939453125e-11
% k 34 = 7.62939453125e-11
% k 35 = 3.814697265625e-11
% k 36 = 3.814697265625e-11
% k 37 = 1.9073486328125e-11
% k 38 = 1.9073486328125e-11
% k 39 = 9.5367431640625e-12
% k 40 = 9.5367431640625e-12
% k 41 = 4.76837158203125e-12
% k 42 = 4.76837158203125e-12
% k 43 = 2.384185791015625e-12
% k 44 = 2.384185791015625e-12
% k 45 = 1.1920928955078126e-12
% k 46 = 1.1920928955078126e-12
% k 47 = 5.960464477539063e-13
% k 48 = 5.960464477539063e-13
% k 49 = 2.9802322387695315e-13
% k 50 = 2.9802322387695315e-13
% k 51 = 1.4901161193847657e-13
% k 52 = 1.4901161193847657e-13
% k 53 = 7.450580596923829e-14

rateconstantswithfrag = [1e-05; 5e-06; 5e-06; 2.5e-06; 2.5e-06; 1.25e-06; 1.25e-06; 6.25e-07; 6.25e-07; 3.125e-07; 3.125e-07; 1.5625e-07;
    1.5625e-07; 7.8125e-08; 7.8125e-08; 3.90625e-08; 3.90625e-08; 1.953125e-08; 1.953125e-08; 9.765625e-09; 9.765625e-09;
    4.8828125e-09; 4.8828125e-09; 2.44140625e-09; 2.44140625e-09; 1.220703125e-09; 1.220703125e-09; 6.103515625e-10;
    6.103515625e-10; 3.0517578125e-10; 3.0517578125e-10; 1.52587890625e-10; 1.52587890625e-10; 7.62939453125e-11; 
    7.62939453125e-11; 3.814697265625e-11; 3.814697265625e-11; 1.9073486328125e-11; 1.9073486328125e-11; 9.5367431640625e-12;
    9.5367431640625e-12; 4.76837158203125e-12; 4.76837158203125e-12; 2.384185791015625e-12; 2.384185791015625e-12;
    1.1920928955078126e-12; 1.1920928955078126e-12; 5.960464477539063e-13; 5.960464477539063e-13; 2.9802322387695315e-13;
    2.9802322387695315e-13; 1.4901161193847657e-13; 1.4901161193847657e-13; 7.450580596923829e-14];

InitialConditions51=[10000,100,100,100,100,100,100,100];%Condition1
InitialConditions52=[5000,100,100,100,100,100,100,100];%Condition2

%in figura 6 cambia solo il fatto che m2=...=m8=1, quindi il resto è uguale
%ma cambiano le condizioni iniziali :
InitialConditions61=[10000,1,1,1,1,1,1,1];%Condition1
InitialConditions62=[5000,1,1,1,1,1,1,1];%Condition2

%in figura 7 cambia solo il fatto che manteniamo le reaction rates constant : 1 per elongation events, 1 per aggregation events, 1 per fragmentation events;
% quindi il resto è uguale e cambiano le condizioni iniziali :
rateconstantsfrag3=[1e-07; 1e-10; 2e-10; 1e-10; 2e-10; 1e-10; 2e-10; 1e-10; 
    2e-10; 1e-10; 2e-10; 1e-10; 2e-10; 1e-10; 2e-10; 1e-10; 2e-10; 1e-10; 
    2e-10; 1e-10; 2e-10; 1e-10; 2e-10; 1e-10; 2e-10; 1e-10; 2e-10; 1e-10; 
    2e-10; 1e-10; 2e-10; 1e-10; 2e-10; 1e-10; 2e-10; 1e-10; 2e-10; 1e-10; 
    2e-10; 1e-10; 2e-10; 1e-10; 2e-10; 1e-10; 2e-10; 1e-10; 2e-10; 1e-10; 
    2e-10; 1e-10; 2e-10; 1e-10; 2e-10; 1e-10];
InitialConditions71=[10000,100,100,100,100,100,100,100];%Condition1
InitialConditions72=[5000,100,100,100,100,100,100,100];%Condition2

%THE INTRODUCTION OF RANDOM SWITCHING DOESN'T REALLY ALTER THE DOMINANCE OF
%CERTAIN REACTION EVENTS : 
%PRIMARY NUCLEATION,
%ELONGATION
%AGGREGATION + MONOMER-DEPENDENT  SECONDARY NUCLEATION

timeEnd=1000;
[timeGrid, avgStates, avgStd, ~, ~] = MonteCarloSimulation(rateconstantswithfrag, Stoichiometry2,InitialConditions51 , timeEnd,100);
[timeGrid2, avgStates2, avgStd2, ~, ~] = MonteCarloSimulation(rateconstantswithfrag, Stoichiometry2,InitialConditions52 , timeEnd,100);
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

    % Assicura che tutti i vettori abbiano la stessa lunghezza
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

    % Assicura che tutti i vettori abbiano la stessa lunghezza
    upperBound = upperBound(:)';
    lowerBound = lowerBound(:)';

    % Ombreggia l'area tra la deviazione standard superiore e inferiore
    fill([timeGridAdjusted, fliplr(timeGridAdjusted)], [upperBound, fliplr(lowerBound)], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    
    title(['Species ' num2str(i)]);
    xlabel('Time');
    ylabel('Population');
    legend('Mean', 'Standard Deviation', 'Location', 'Best');
end

% Aggiungo i risultati della seconda simulazione alle figure esistenti
figure(1); % Seleziona la prima figura
for i = 1:4
    subplot(2, 2, i);
    hold on;
    % Plotta la dinamica media della seconda simulazione
    plot(timeGrid2, squeeze(avgStates2(1, i, :)), 'g', 'LineWidth', 2); % Usa il colore verde per distinguere la seconda simulazione
    
    % Calcola i limiti superiori e inferiori dell'area ombreggiata per la seconda simulazione
    upperBound2 = squeeze(avgStates2(1, i, :) + avgStd2(1, i, :));
    lowerBound2 = squeeze(avgStates2(1, i, :) - avgStd2(1, i, :));
    % Regola la lunghezza dei vettori per evitare errori
    timeGridAdjusted = linspace(0, timeEnd, length(upperBound));

    % Assicura che tutti i vettori abbiano la stessa lunghezza
    upperBound2 = upperBound2(:)';
    lowerBound2 = lowerBound2(:)';
    % Ombreggia l'area tra la deviazione standard superiore e inferiore per la seconda simulazione
    fill([timeGridAdjusted, fliplr(timeGridAdjusted)], [upperBound2, fliplr(lowerBound2)], 'm', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    legend('Mean (Condition 1)', 'Std Dev (Condition 1)', 'Mean (Condition 2)', 'Std Dev (Condition 2)', 'Location', 'Best');

    
end

figure(2); % Seleziona la seconda figura
for i = 5:size(avgStates2, 2)
    subplot(2, 2, i-4);
    hold on;
    % Plotta la dinamica media della seconda simulazione
    plot(timeGrid2, squeeze(avgStates2(1, i, :)), 'g', 'LineWidth', 2); % Usa il colore verde per distinguere la seconda simulazione
    
    % Calcola i limiti superiori e inferiori dell'area ombreggiata per la seconda simulazione
    upperBound2 = squeeze(avgStates2(1, i, :) + avgStd2(1, i, :));
    lowerBound2 = squeeze(avgStates2(1, i, :) - avgStd2(1, i, :));
    % Regola la lunghezza dei vettori per evitare errori
    timeGridAdjusted = linspace(0, timeEnd, length(upperBound));

    % Assicura che tutti i vettori abbiano la stessa lunghezza
    upperBound2 = upperBound2(:)';
    lowerBound2 = lowerBound2(:)';
    % Ombreggia l'area tra la deviazione standard superiore e inferiore per la seconda simulazione
    fill([timeGridAdjusted, fliplr(timeGridAdjusted)], [upperBound2, fliplr(lowerBound2)], 'm', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    legend('Mean (Condition 1)', 'Std Dev (Condition 1)', 'Mean (Condition 2)', 'Std Dev (Condition 2)', 'Location', 'Best');
end

% Inizializza i vettori per il momento zero e il primo momento
zerothMoment = zeros(1, size(avgStates, 3));
firstMoment = zeros(1, size(avgStates, 3));

% Calcola il momento zero e il primo momento per ogni istante di tempo
for t = 1:size(avgStates, 3)
    zerothMoment(t) = sum(avgStates(1, 2:end, t));
    for s = 2:size(avgStates, 2) % Itera da 2 a size(avgStates, 2) per escludere la prima specie
        firstMoment(t) = firstMoment(t) + sum(avgStates(1, s, t) * s);
    end
end

% Crea una nuova figura
figure(3);

% Plotta il momento zero
subplot(3, 1, 1);
hold on ;
plot(1:size(avgStates, 3), zerothMoment, 'b', 'LineWidth', 2);
title('Zeroth Moment');
xlabel('Time');
ylabel('Total Population (excluding species 1)');

% Plotta il primo momento
subplot(3, 1, 2);
plot(1:size(avgStates, 3), firstMoment, 'r', 'LineWidth', 2);
title('First Moment');
xlabel('Time');
ylabel('First Moment (Mass Concentration) (excluding species 1)');

%Plotta il rapporto M/T
subplot (3,1,3)
plot(1:size(avgStates, 3), firstMoment./zerothMoment, 'g', 'LineWidth', 2);
title('First Moment/Zeroth Moment');
xlabel('Time');
