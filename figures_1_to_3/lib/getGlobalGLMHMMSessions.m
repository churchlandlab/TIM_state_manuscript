function sessiondates = getGlobalGLMHMMSessions(glmPath)
%dates = cellstr(load(glmPath).model_training_sessions);
load(glmPath);
dates = cellstr(model_training_sessions);
%session_animals = cellstr(load(glmPath).mouse);
session_animals = cellstr(mouse);
animals = unique(session_animals);
for i = 1:length(animals)
    sessiondates{i}  = dates(ismember(session_animals,animals(i)));
end