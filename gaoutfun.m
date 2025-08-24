
%% function file: save below code in separate file and name it as "gaoutfun"

function [state,options,optchanged] = gaoutfun(options,state,flag)
persistent  history scores
optchanged = false;
switch flag
    case 'init'
       
        history(:,:,1) = state.Population;
        scores(:,1) = state.Score;
        assignin('base','gapopulationhistory',history);
        assignin('base','allscore',scores);
    case 'iter'
        % Update the history every 1 generations.
        if rem(state.Generation,1) == 0
            ss = size(history,3);
            history(:,:,ss+1) = state.Population;
            scores(:,ss+1) = state.Score;
            assignin('base','gapopulationhistory',history);
            assignin('base','allscore',scores);
        end
    case 'done'
        % Include the final population in the history.
        ss = size(history,3);
        history(:,:,ss+1) = state.Population;
        scores(:,ss+1) = state.Score;
        assignin('base','gapopulationhistory',history);
        assignin('base','allscore',scores);
end