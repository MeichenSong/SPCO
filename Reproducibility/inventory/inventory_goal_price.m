% function gp=inventory_goal_price(alpha, salvaged_value, cost, price)
% % function of g(p) for multiplication model with dp=30-price, price in [cost,100]
% dp=20-price;
% r=(1-alpha)*(price-cost)/(price-salvaged_value);
% funs=@(x) x.^2.*(1-x/2)/(4*beta(2,2));
% up=2*betainv(r,2,2);
% gp=(price-salvaged_value)*dp*integral(funs,0,2*betainv(r,2,2))/(1-alpha);
% gp=((price-salvaged_value)*dp/(1-alpha))*((up^3)/3-(up^4)/8)/(4*beta(2,2));

%%%%%%%%%% CVaR %%%%%%%%
function gp=inventory_goal_price(alpha, theta)
% function of g(p) for multiplication model with dp=30-price, price in [cost,100]
b=(theta(2)-1)*(100-theta(2)); a = -3*theta(1);
if alpha==0.05
    quantile_ori=0.86465;
else
    quantile_ori=0.941097;
end
upper_int = 2*b*quantile_ori+a;
funs=@(x) (x/(4*b*beta(2,2))).*((x-a)/b).*(1-(x-a)/(2*b));
gp=integral(funs,a,upper_int)/(1-alpha);