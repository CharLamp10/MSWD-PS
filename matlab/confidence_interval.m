function [mean_value,lower_bound,upper_bound,margin_of_error] = confidence_interval(data,varargin)

if isempty(varargin)
    confidence_level = 0.95;
else
    confidence_level = varargin{1};
end
mean_value = mean(data,'omitnan');

std_error = std(data,'omitnan')./sqrt(repmat(size(data,1),[1,size(data,2)]));

degrees_of_freedom = repmat(size(data,1) - 1,[1,size(data,2)]);
alpha = 1 - confidence_level;
t_critical = tinv(1 - alpha / 2, degrees_of_freedom);

margin_of_error = t_critical .* std_error;

% Calculate the 95% confidence interval
lower_bound = mean_value - margin_of_error;
upper_bound = mean_value + margin_of_error;


end