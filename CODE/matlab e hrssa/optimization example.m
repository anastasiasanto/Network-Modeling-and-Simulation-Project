%calcolo della distanza euclidea tra il primo momento sperimentale e
%simulato

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
rateconstants = [1e-05; 5e-06; 2.5e-06; 1.25e-06; 6.25e-07;3.125e-07; 1.5625e-07; 7.8125e-08; 3.90625e-08; 1.953125e-08;
    9.765625e-09; 4.8828125e-09; 2.44140625e-09; 1.220703125e-09; 6.103515625e-10; 3.0517578125e-10; 1.52587890625e-10; 
    7.62939453125e-11; 3.814697265625e-11; 1.9073486328125e-11; 9.5367431640625e-12; 4.76837158203125e-12; 2.384185791015625e-12;
    1.1920928955078126e-12; 5.960464477539063e-13; 2.9802322387695315e-13; 1.4901161193847657e-13];

InitialConditions31=[10000,1,1,1,1,1,1,1];%Condition1
InitialConditions32=[5000,1,1,1,1,1,1,1];%Condition2
% Dati sperimentali (da sostituire con i tuoi dati reali)
    experimental_moments = [0.14,0.19,0.23,0.27,0.32,0.36,0.39,0.42,0.45,0.48,0.49,0.51,0.52,0.56,0.55,
        0.56,0.59,0.61,0.62,0.64,0.67,0.67,0.72,0.74,0.71,0.75,0.78,0.79,0.81,0.8,0.82,0.81,0.87,0.88,0.9,0.88];
    texp = [76.19,128.91,209.94,299.07,388.21,477.34,566.45,655.57,744.67,833.77,922.86,1018.88,1101.01,1185.07,
        1279.18,1368.24,1457.35,1546.43,1635.52,1720.56,1805.62,1886.57,1955.48,2032.42,2121.42,2190.31,2259.18,2348.26,
        2437.34,2526.39,2619.83,2700.46,2776.46,2843.88,2931.35,2992.03] 
    
    function distance = calculate_distance(rate_constants, stoichiometrymatrix, initialconditions)
    % Simula il modello utilizzando le rate constants fornite
    timeEnd=1000;
    [timeGrid, avgStates, avgStd, ~, ~] = MonteCarloSimulation(rateconstants, Stoichiometrynofrag,InitialConditions31 , timeEnd,100);
    for t = 1:size(avgStates, 3)
        for s = 2:size(avgStates, 2) % Itera da 2 a size(avgStates, 2) per escludere la prima specie
            firstMoment(t) = firstMoment(t) + sum(avgStates(1, s, t) * s); %easier formula
        end
    end
    % Calcola la distanza euclidea tra i momenti simulati e quelli sperimentali
    distance = norm(simulated_moments - experimental_moments);
end

distance = calculate_distance(initial_rate_constants);

%The computation of the moment should be done with the proper formula given
%in the paper, but in this context we were not able to apply the
%mathematical computations in order ot write it and extract M(t9 as a
%function of k+, k-, kn, koff, k2.
%distance = norm(simulatedfirstmoment - experimentafirstmoment) --> we can
%express this as a function of k+,kn in the case of aggregation only, and
%k- ,k2 and koff in the case in which fragmentation is allowed. We can put this
%function into an optimization algorithm and use it for the fitness
% evaluation such that for different values of k ( generated initially at random)
% we can calculate the euclidean distance between experimental and simulated data.
% The purpose of the optimization is to minimize this distance.
% In particular since it is a discrete function we have to 
% use Genetic algorithm or genetic programming , in which we have the
%following parameters to set : FROM PYTHON CODE 

args["4"] = 10 # Number of dimensions of the search space, so number of variables we want to optimize
args["gaussian_stdev"] = 0.5 # Standard deviation of the Gaussian mutations
args["mutation_rate"] = 0.5 # fraction of loci to perform mutation on
args["tournament_size"] = 2 
args["num_elites"] = 1 # number of elite individuals to maintain in each gen
args["pop_size"] = 20 # population size
args["pop_init_range"] = [0.0000000001, 0.0001] # Range for the initial population
args["max_generations"] = 50 # Number of generations of the GA

%Using the package 'inspyred' we can be able to optimize these values.
