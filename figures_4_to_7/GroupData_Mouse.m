% Divide trials into different plotting groups based on divide_mode
function [InputMatrix] = GroupData_Mouse(InputMatrix, divide_mode, raw_data)

InputMatrix(:, end+1) = 0;   % Mark for dividing trials

%% 'stimulus'
if isequal(divide_mode, 'stimulus') == 1
    for a = 1 : length(raw_data.CorrectSide)
        InputMatrix(a, end) = raw_data.CorrectSide(a); 
    end


%% 'opto'
elseif isequal(divide_mode, 'opto') == 1
    if isfield(raw_data, 'optoType') == 1
        for a = 1 : length(raw_data.optoType)
            if ~isnan(raw_data.optoType(a))
                InputMatrix(a, end) = 1;
            else
                InputMatrix(a, end) = 0;
            end
        end
        
    else
        InputMatrix(:, end) = 0;
    end


%% 'optotype'
elseif isequal(divide_mode, 'optotype') == 1
    if isfield(raw_data, 'optoType') == 1
        for a = 1 : length(raw_data.optoType)
            if ~isnan(raw_data.optoType(a))
                InputMatrix(a, end) = raw_data.optoType(a);
            else
                InputMatrix(a, end) = 0;    % No opto stimulation
            end
        end
        
    else
        InputMatrix(:, end) = 0;
    end
    
    
%% 'output'
elseif isequal(divide_mode, 'outcome') == 1
    for a = 1 : length(raw_data.Rewarded)
        if raw_data.DidNotChoose(a) == 1
            InputMatrix(a, end) = 0;	% No response
        elseif raw_data.DidNotChoose(a) == 0 && raw_data.Rewarded(a) == 1
            InputMatrix(a, end) = 1;    % Correct
        elseif raw_data.DidNotChoose(a) == 0 && raw_data.Rewarded(a) == 0
            InputMatrix(a, end) = 2;    % Wrong
        end
    end
    
    
%% 'responseside'
elseif isequal(divide_mode, 'responseside') == 1
    for a = 1 : length(raw_data.ResponseSide)
        if isnan(raw_data.ResponseSide(a))
            InputMatrix(a, end) = 0;	% No response
        elseif ~isnan(raw_data.ResponseSide(a)) && raw_data.ResponseSide(a) == 1
            InputMatrix(a, end) = 1;    
        elseif ~isnan(raw_data.ResponseSide(a)) && raw_data.ResponseSide(a) == 2
            InputMatrix(a, end) = 2;    
        end
    end
    
    
%% 'formeroutcome'
elseif isequal(divide_mode, 'formeroutcome') == 1
    InputMatrix(1, end) = 0;    % this 1st trial doesn't have a former trial. Setting is as 'no response'
    for a = 2 : length(raw_data.Rewarded)
%         if raw_data.DidNotChoose(a-1) == 1
%             InputMatrix(a, end) = 0;	% No response
        if raw_data.Rewarded(a-1) == 1
            InputMatrix(a, end) = 2;    % Correct
        elseif raw_data.Rewarded(a-1) == 0
            InputMatrix(a, end) = 1;    % Wrong/No response
        end
    end
    
    
%% 'formerresponseside'
elseif isequal(divide_mode, 'formerresponseside') == 1
    InputMatrix(1, end) = 0;    % this 1st trial doesn't have a former trial. Setting is as 'no response'
    for a = 2 : length(raw_data.ResponseSide)
        if isnan(raw_data.ResponseSide(a-1))
            InputMatrix(a, end) = 0;	% No response
        elseif ~isnan(raw_data.ResponseSide(a-1)) && raw_data.ResponseSide(a-1) == 1
            InputMatrix(a, end) = 1;
        elseif ~isnan(raw_data.ResponseSide(a-1)) && raw_data.ResponseSide(a-1) == 2
            InputMatrix(a, end) = 2;
        end
    end

    
%% 'random'
elseif isequal(divide_mode, 'random') == 1
    for a = 1 : size(InputMatrix, 1)
        if rand < 0.5
            InputMatrix(a,end) = 0;
        else
            InputMatrix(a,end) = 1;
        end
    end

    
%% 'optotype_classifier'    This is especially designed for function "LinearClassifierPlus"
elseif isequal(divide_mode, 'optotype_classifier') == 1
    for a = 1 : length(raw_data.optoType)
        if ~isnan(raw_data.optoType(a))
            InputMatrix(a, end) = raw_data.optoType(a)*raw_data.optoArea(a);      
        else
            InputMatrix(a, end) = raw_data.optoArea(a);    % No opto stimulation
        end
    end
    
end
    
end


