% function [price_opt,order_opt]=optimal_price_order(alpha,cost,salvaged_value)
% lb=cost;ub=20;
% 
% 
% A = [];
% b = [];
% Aeq = [];
% beq = [];
% options = optimoptions(@ga,'PlotFcn',{@gaplotbestf,@gaplotstopping},'FunctionTolerance',1e-6);
% % options = optimoptions(@ga,'PlotFcn',{@gaplotbestf,@gaplotstopping},'MaxStallGenerations',1000,'FunctionTolerance',1e-6,'CrossoverFraction',0.8,'Display','off');
% % [theta_opt,CVaR_opt] = fmincon(@(opt_variable) -RockGoal_Low(input_samples, opt_variable, quantile_level), theta0,A,b,Aeq,beq,lb,ub);
% price_opt = ga(@(price) -inventory_goal_price(alpha, salvaged_value, cost, price), 1,A,b,Aeq,beq,lb,ub,[],options);
% % price_opt = fmincon(@(price) -inventory_goal_price(alpha, salvaged_value, cost, price), cost,A,b,Aeq,beq,lb,ub);
% % price_opt = patternsearch(@(price) -inventory_goal_price(alpha, salvaged_value, cost, price), 70,A,b,Aeq,beq,lb,ub);
% r=(1-alpha)*(price_opt-cost)/(price_opt-salvaged_value);
% order_opt = (20-price_opt)*2*betainv(r,2,2);
% % i=1;
% % price_ca=[4:0.01:100];
% % for i=1:size(price_ca,2)
% %     price=price_ca(i);
% %     goal(i)=inventory_goal_price(alpha, salvaged_value, cost, price);
% % %     i=i+1;
% % end
% % [a,b]=max(goal);
% % price_opt=price_ca(b);
% % r=(1-alpha)*(price_opt-cost)/(price_opt-salvaged_value);
% % order_opt = (100-price_opt)*2*betainv(r,2,2);


%%%%%%%%%%%% opt for CVaR %%%%%%%5
% 
% function theta_opt=optimal_price_order(alpha)
% lb=[0,4];ub=[100,100];
% theta0=[20,20];
% 
% A = [-1,-2];
% b = -200;
% Aeq = [];
% beq = [];
% options = optimoptions(@ga,'FunctionTolerance',1e-6);
% % theta_opt = ga(@(theta) -inventory_goal_price(alpha, theta), 2,A,b,Aeq,beq,lb,ub,[],options);
% theta_opt = fmincon(@(theta) -inventory_goal_price(alpha, theta), theta0,A,b,Aeq,beq,lb,ub);
% % theta_opt = patternsearch(@(theta) -inventory_goal_price(alpha, theta), theta0,A,b,Aeq,beq,lb,ub);
alpha=0.05;
theta2=[50:0.01:100]';
theta1=200-2*theta2;
theta_a(:,1)=theta1;
theta_a(:,2)=theta2;
for i=1:size(theta2,1)
    theta=theta_a(i,:);
    gp(i)=inventory_goal_price(alpha, theta);
end
plot([1:1:size(theta2,1)],gp)
[value,index]=max(gp);