function [timeGrid, avgStates, avgStd, allReactionsOccurred, allReactionTimes] = MonteCarloSimulation(rateConstants, stoichiometry, initialState, timeEnd, numSimulations)
    numSpecies = length(initialState);
    allReactionsOccurred = [];
    allReactionTimes = [];
    timeGrid = linspace(0, timeEnd, 1000);
    statesInterpolated = zeros(numSimulations, numSpecies, length(timeGrid));

    for i = 1:numSimulations
        [times, states, reactionsOccurred, reactionTimes] = RSSA(rateConstants, stoichiometry, initialState, timeEnd);

        for j = 1:numSpecies
            statesInterpFunc = interpolate(times, states(j, :), timeGrid);
            statesInterpolated(i, j, :) = statesInterpFunc(timeGrid);
        end

        allReactionsOccurred{i} = reactionsOccurred;
        allReactionTimes{i} = reactionTimes;
    end

    avgStates = mean(statesInterpolated, 1);
    avgStd = std(statesInterpolated, 1);
end
function interpFunc = interpolate(x, y, xq)
    if length(x) > 1
        interpFunc = griddedInterpolant(x, y, 'linear', 'nearest');
        interpFunc = @(xq) interpFunc(xq);
    else
        % If x is a single value, use linear interpolation directly
        interpFunc = @(xq) interp1(x, y, xq, 'linear', 'nearest');
    end
end
% Require: rate constants (list): The rate constants for the reactions.
% Require: stoichiometry (list): The stoichiometry of the reactions. A 2D
% list where the i-th list contains the stoichiometric coefficients of the i-th
% reaction.
% Require: initial state (list): The initial numbers of molecules for each
% species.
% Require: time end (float): The simulation time end.
% Require: num simulations (int): The number of simulations to run.
