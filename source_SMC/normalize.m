% normalize a log of weight 
function weight = normalize(log_weight)
log_weight = reshape(log_weight, [numel(log_weight), 1]); 
 weight = exp(log_weight' - max(log_weight'));
 weight = weight ./ sum(weight);
end