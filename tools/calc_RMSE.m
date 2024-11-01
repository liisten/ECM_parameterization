function rmse = calc_RMSE(a,b)

rmse = sqrt(mean((a(:)-b(:)).^2));

end