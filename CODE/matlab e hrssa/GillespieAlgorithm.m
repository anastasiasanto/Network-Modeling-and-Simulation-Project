function [times, states, reactionsOccurred, reactionTimes] = GillespieAlgorithm(rateConstants, stoichiometry, initialState, timeEnd)
    numReactions = length(rateConstants);
    numSpecies = length(initialState);
    state = initialState;
    time = 0;
    times = [time];
    states = [state];
    reactionsOccurred = [];
    reactionTimes = [];

    while time < timeEnd
        propensities = zeros(1, numReactions);

        for i = 1:numReactions
            propensity = rateConstants(i);
            for j = 1:numSpecies
                if stoichiometry(i, j) < 0
                    propensity = propensity * state(j) ^ abs(stoichiometry(i, j));
                end
            end
            propensities(i) = propensity;
        end

        totalPropensity = sum(propensities);

        if totalPropensity <= 0
            break;
        end

        time = time + exprnd(1/totalPropensity);
        reactionIndex = randsample(1:numReactions, 1, true, propensities);
        
        state = state + stoichiometry(reactionIndex, :);
        times = [times, time];
        states = [states; state];
        reactionsOccurred = [reactionsOccurred, reactionIndex];
        reactionTimes = [reactionTimes, time];
    end

    % Transpose the states matrix for compatibility with the pseudocode
    states = states';
end
