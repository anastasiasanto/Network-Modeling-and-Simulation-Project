function [times, states, reactionsOccurred, reactionTimes] = NextReactionMethod(rateConstants, stoichiometry, initialState, timeEnd)
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
        randV = rand(1,initialLenght);
        nextReactionTime = zeros(size(rateConstants));
        % Next Reaction Method: select the next reaction and time 
        for i = 1:numReactions 
            nextReactionTime(i) = propensity(i)*log(1/randV(j));
        end
        mint=min(nextReactionTime);
        reactionIndex=find(nextReactionTime==mint);
        % Definizione di un vettore di esempio

        % Update state
        state = state + stoichiometry(reactionIndex, :);

        % Update time
        time = time + nextReactionTime;

        % Store results
        times = [times, time];
        states = [states; state];
        reactionsOccurred = [reactionsOccurred, reactionIndex];
        reactionTimes = [reactionTimes, time];
    end

    % Transpose the states matrix for compatibility with the pseudocode
    states = states';
end
