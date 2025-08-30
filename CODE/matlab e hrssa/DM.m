function [times, states, reactionsOccurred, reactionTimes] = DM(rateConstants, stoichiometry, initialState, timeEnd)
    numReactions = length(rateConstants);
    numSpecies = length(initialState);
    state = initialState;
    time = 0;
    times = [time];
    states = [state];
    reactionsOccurred = [];
    reactionTimes = [];

    while time < timeEnd
        % Calculate reaction propensities
        propensities = zeros(1, numReactions);
        for i = 1:numReactions
            propensities(i) = rateConstants(i);
            for j = 1:numSpecies
                if stoichiometry(i, j) < 0
                    propensities(i) = propensities(i) * state(j) ^ abs(stoichiometry(i, j));
                end
            end
        end

        % Calculate total propensity
        totalPropensity = sum(propensities);

        % Generate a candidate reaction time
        candidateReactionTime = exprnd(1/totalPropensity);

        % Generate a random number to decide which reaction occurs
        reactionSelector = rand() * totalPropensity;

        % Initialize variables for candidate reaction
        candidateReactionIndex = 0;
        cumulativePropensity = 0;

        % Find the candidate reaction
        for i = 1:numReactions
            cumulativePropensity = cumulativePropensity + propensities(i);
            if cumulativePropensity >= reactionSelector
                candidateReactionIndex = i;
                break;
            end
        end

        % Check if the candidate reaction time is within the simulation time
        if (time + candidateReactionTime) > timeEnd
            break;
        end

        % Update the system state based on the candidate reaction
        state = state + stoichiometry(candidateReactionIndex, :);

        % Update time
        time = time + candidateReactionTime;

        % Store results
        times = [times, time];
        states = [states; state];
        reactionsOccurred = [reactionsOccurred, candidateReactionIndex];
        reactionTimes = [reactionTimes, time];
    end

    % Transpose the states matrix for compatibility with the pseudocode
    states = states';
end
