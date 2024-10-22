function RMSE = rmse(output, target)

RMSE = sqrt(sum(abs(output - target).^2)/length(output));

end