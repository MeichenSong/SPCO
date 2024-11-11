%%%% inventory with unfixed price, alpha=0.05/0.01
function [theta,q,CVaR,CVaR_star_hat]=SPCO_inventory_one(K,alpha, WarmUp, bEnlarge, cEnlarge,skiprow) 
q=zeros(K,1);
cost=4; salvaged_value=1; %price=10;
% d=2;
% initialize theta
theta=[10,10];
b_upper=[100,100];
b_lower=[0,cost];
C=ceil(K*WarmUp);
batch_size=10000;

for k=2:K
   

%     bk=bEnlarge*(1)/(k-1);
% %     ck=cEnlarge*(1)/(k-1)^(0.167);
%     ck=cEnlarge*(1)/((max(k-0.1*K,1))^(0.167));
%     bk=bEnlarge*(1)/(k+10*bEnlarge);
%     ck=cEnlarge*(1)/(k+10*cEnlarge)^(0.167);
%     gammak=1*C/(k-1)^(5/9);
%     
    %%%%%%
    bk=bEnlarge*(1)/((k+10*bEnlarge)^0.99);
    ck=cEnlarge*(1)/((k+10*cEnlarge)^0.1429);
    gammak=C/(k-1)^(0.5714);
    
    delta=2*(rand(1,2)>0.5)-1;
    
    thetan=theta(k-1,:);
    thetap=thetan+ck*delta; % pertubation size for price and order; magnitute
    thetam=thetan-ck*delta;
    
    CVaR_hat(k-1)=inventory_goal_price(alpha, thetan);
    %%% approximate CVaR(thetak) %%%%%%
%     X_sample=2*betarnd(2,2,1,batch_size);
%     demand_sample=(100-thetan(2))*X_sample;
%     Y_sample=(thetan(2)-cost)*thetan(1)-(thetan(2)-salvaged_value)*max(thetan(1)-demand_sample,0);
%     Y_sample_sort = sort(Y_sample);
%     q_hat = Y_sample_sort(ceil((1-alpha)*batch_size));
%     CVaR_hat(k-1) = sum(Y_sample.*(Y_sample<=q_hat))/((1-alpha)*batch_size);
%     clear X_sample demand_sample Y_sample Y_sample_sort q_hat
    %%%%%%%%%%%

    % sample input
%     X=gamrnd(2,2,1,3); % 1*3
    X=2*betarnd(2,2,1,3);
    demandn=(100-thetan(2))*X(1);
    demandp=(100-thetap(2))*X(2);
    demandm=(100-thetam(2))*X(3);
    Y=(thetan(2)-cost)*thetan(1)-(thetan(2)-salvaged_value)*max(thetan(1)-demandn,0); %profit
    Yp=(thetap(2)-cost)*thetap(1)-(thetap(2)-salvaged_value)*max(thetap(1)-demandp,0);
    Ym=(thetam(2)-cost)*thetam(1)-(thetam(2)-salvaged_value)*max(thetam(1)-demandm,0);
    
    thetatmp = theta(k-1,:)+ bk*((min(q(k-1),Yp)-min(q(k-1),Ym))./(2*(1-alpha)*ck*delta));
    q(k)=q(k-1)+gammak*(1-alpha-(Y<=q(k-1))); 
    
    %projection
%     k
    if thetatmp(1)-0.5*thetatmp(2)>=75 && thetatmp(2)<=50
%         1
        theta(k,:)=[100,50];
    elseif thetatmp(1)>100 && 50<thetatmp(2) && thetatmp(2)<=100
%         2
        theta(k,:)=[100,thetatmp(2)]; 
    elseif thetatmp(1)>100 && thetatmp(2)>100
%         3
        theta(k,:)=[100,100];
    elseif thetatmp(1)>0 && thetatmp(1)<=100 && thetatmp(2)>100
%         4
        theta(k,:)=[thetatmp(1),100];
    elseif thetatmp(1)-0.5*thetatmp(2)<=-50 && thetatmp(1)<=0
%         5
        theta(k,:)=[0,100];
    elseif thetatmp(1)-0.5*thetatmp(2)>-50 && thetatmp(1)-0.5*thetatmp(2)<75 && thetatmp(1)+2*thetatmp(2)<200 
%         6
        theta(k,1)=40-0.4*(thetatmp(2)-2*thetatmp(1));
        theta(k,2)=80+0.2*(thetatmp(2)-2*thetatmp(1));
    else
        theta(k,:)=thetatmp;
    end
        
    
%     thetatmp = thetatmp.*(thetatmp<=b_upper)+b_upper.*(thetatmp>b_upper);
%     theta(k,:) = thetatmp.*(thetatmp>=b_lower)+b_lower.*(thetatmp<b_lower);
end
%%% approximate CVaR_star %%%%%%
switch alpha
case 0.05
    theta_star(1) = 92.7288;
    theta_star(2) = 53.6356;
case 0.01
    theta_star(1) = 92.9436;
    theta_star(2) = 53.5282;
end
% X_sample=2*betarnd(2,2,1,batch_size);
% demand_sample=(100-theta_star(2))*X_sample;
% Y_sample=(theta_star(2)-cost)*theta_star(1)-(theta_star(2)-salvaged_value)*max(theta_star(1)-demand_sample,0);
% Y_sample_sort = sort(Y_sample);
% q_hat = Y_sample_sort(ceil((1-alpha)*batch_size));
% CVaR_star_hat = sum(Y_sample.*(Y_sample<=q_hat))/((1-alpha)*batch_size);
CVaR_star_hat =inventory_goal_price(alpha, theta_star);
%%%%%%%%%%%
q=q(1:skiprow:end);
theta=theta(1:skiprow:end,:);
CVaR = CVaR_hat(1:skiprow:end);
% subplot(3,1,1)
% plot([1:skiprow:K],theta(:,1),'k-', 'LineWidth',1.5)
% hold on
% plot([1:skiprow:K],theta_star(1)*ones(size([1:skiprow:K])),'r-', 'LineWidth',1.5)
% plot([1:skiprow:K],theta(:,2),'b-', 'LineWidth',1.5)
% plot([1:skiprow:K],theta_star(2)*ones(size([1:skiprow:K])),'r-', 'LineWidth',1.5)
% legend("order","price")
% ylabel('$\theta_k$','interpreter','latex','FontSize',15);
% title('$\alpha=0.95$, Fixed Price','interpreter','latex','FontSize',15);
% subplot(3,1,2)
% plot([1:skiprow:K],q,'k-', 'LineWidth',1.5)
% ylabel('$q_k$','interpreter','latex','FontSize',15);
% subplot(3,1,3)
% plot([1:skiprow:K],CVaR,'b-', 'LineWidth',1.5)
% hold on
% plot([1:skiprow:K],CVaR_star_hat*ones(size(CVaR)),'r-', 'LineWidth',0.5)

    
    
    
    
       
 
 
 