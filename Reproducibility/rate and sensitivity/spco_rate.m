% this file is created on May 27 2024:
% 1. add an exp for checking the optimal rate

function [CVaR,CVaR_star,dist_spco,theta]=spco_rate(K,alpha,experi, WarmUp, bEnlarge, cEnlarge,skiprow,beta,gamma,c) 
    d=20;
    nk_2=1;
    CVaR=zeros(K,1);q=zeros(K,1);
    dist_spco=zeros(K,1);
    
    % initialize theta
    switch experi
    case 2 
        theta = [1:d]'.*rand(d,K);
        theta_star = 0.5*[1:d].';
        theta_lower=zeros(d,1);
        theta_upper=[1:d].';
  
    case 1
        theta = [1:d]'+rand(d,K)-0.5;
        theta_star = [1:d].';
        theta_upper = [1:d]'+0.5;
        theta_lower=[1:d]'-0.5;
    end
    
    C=ceil(K*WarmUp);
    for k=2:K
%         N=N+nk_1+2*nk_2;
       
%         bk=bEnlarge*(2*C)/(k-1+C);
%         ck=cEnlarge*((2*C)^(0.167))/(k-1+C)^(0.167);
%         bk=bEnlarge*(1)/(k-1);
%%%%%%% origin parameter order
%         bk=bEnlarge*(1)/((k+10*bEnlarge)^0.99);
%         ck=cEnlarge*(1)/((k+10*cEnlarge)^0.167);
%         gammak=C/(k-1)^(5/9);

%%%%%%%%%%%%%%%% update order
%         bk=bEnlarge*(1)/((k)^0.99);
%         ck=cEnlarge*(1)/((k)^0.1429);
        
%         if opt_rate_index ==1
% %             beta = 0.99;
%         bk=bEnlarge*(1)/((k)^beta);
%         ck=cEnlarge*(1)/((k)^0.1429);
% %             bk=bEnlarge*(1)/((k+10*bEnlarge)^0.99);
% %             ck=cEnlarge*(1)/((k+10*cEnlarge)^0.1429);
%             gammak=C/(k-1)^(0.5714);
%             delta=2*(rand(d,1)>0.5)-1;
%         else
%             beta = 0.8;
%             bk=bEnlarge*(1)/((k)^beta);
%             ck=cEnlarge*(1)/((k)^0.1429);
%             gammak=C/(k-1)^(0.5714);
%             delta=2*(rand(d,1)>0.5)-1;
%         end
        bk=bEnlarge*(1)/((k)^beta);
        ck=cEnlarge*(1)/((k)^c);
        gammak=C/(k)^(gamma);
        delta=2*(rand(d,1)>0.5)-1;

%         if experi == 4
%             delta=2*pi*(rand(d,1)>0.5)/d-pi/d;
%         end
        
        thetan=theta(:,k-1);
        thetap=thetan+ck*delta;
        thetam=thetan-ck*delta;
        
        % sample input
        p = rand(nk_2,3);
        X=sqrt(2) * erfinv(2 * p - 1);

        switch experi
        case 2
            Y=X(1,1)+h2_rate(thetan,d); Yp=X(:,2)+h2_rate(thetap,d); Ym=X(:,3)+h2_rate(thetam,d); 

        case 1
            Y=X(1,1)*h1_rate(thetan,d); Yp=X(:,2)*h1_rate(thetap,d); Ym=X(:,3)*h1_rate(thetam,d); 
        end 
      
        q(k)=q(k-1)+gammak*(1-alpha-(Y<=q(k-1)));
        
        thetatmp = theta(:,k-1)-bk*(mean((max(q(k-1),Yp))-mean(max(q(k-1),Ym)))./(2*alpha*ck*delta));
        
%         indi = -(1-alpha-(Y<=q(k-1))<=0)+ (1-alpha-(Y<=q(k-1))>=0);
%         q(k)=q(k-1)+1*sign(1-alpha-(Y<=q(k-1)));
        
        thetatmp= thetatmp.*(thetatmp<=theta_upper)+theta_upper.*(thetatmp>theta_upper);
        theta(:,k)= thetatmp.*(thetatmp>=theta_lower)+theta_lower.*(thetatmp<theta_lower);
        switch experi
        case 1
            CVaR(k)=2.0627*h1_rate(theta(:,k),d);
        case 2
            CVaR(k)=2.0627+h2_rate(theta(:,k),d);
        end
        dist_spco(k) = vecnorm(theta(:,k)-theta_star);
    end

    switch experi
        case 1
            CVaR_star = h1_rate(theta_star,d) *2.0627;
        case 2
            CVaR_star = h2_rate(theta_star,d) +  2.0627;
    end

    
    dist_spco=dist_spco(1:skiprow:end);  
    log_dist=log(dist_spco);        
    CVaR = CVaR(1:skiprow:end);
    theta = theta(:,1:skiprow:end);
%     k = [1:skiprow:K]; 
%     Logk = log(k);
%     x=linspace(0,max(Logk),20);
%     y_SPCO=interp1(Logk,log_dist,x);
%     q = q(1:skiprow:end);
%     theta1 = theta(1,1:skiprow:end);
%     theta2 = theta(2,1:skiprow:end);
%     theta = theta(:,1:skiprow:end);
%     regress_begin =  ceil(0.7*size(x,2));
%     regress_end = ceil(0.95*size(x,2));
%     X = [ones(size(x));x];
%     Beta_SPCO = regress(y_SPCO(regress_begin:end).',X(:,regress_begin:end).');
%     
%     theta3 = theta(3,1:skiprow:end);
%     theta4 = theta(4,1:skiprow:end);
%     theta5 = theta(5,1:skiprow:end);
%     theta6 = theta(6,1:skiprow:end);
%     theta7 = theta(7,1:skiprow:end);
%     theta8 = theta(8,1:skiprow:end);
%     theta9 = theta(9,1:skiprow:end);
%     theta10 = theta(10,1:skiprow:end);
%     subplot(4,1,4);
%     plot(Logk,log_dist,'bo')
%     hold on
%     plot([0:max(Logk)+1],Beta_SPCO(2)*[0:max(Logk)+1]+Beta_SPCO(1),'k-','LineWidth',1.5);
%     Rlinename_opt = sprintf('y=%.2fx + %.2f', Beta_SPCO(2), Beta_SPCO(1));
%     legend(Rlinename_opt, 'interpreter','latex','Location','best');
%     subplot(4,1,1);
%     plot(3*[1:skiprow:K],CVaR,'k-', 'LineWidth',1.5);
%     ylabel('$CVaR(\theta_k)$','interpreter','latex','FontSize',15);
%     hold on
%     plot(3*[1:skiprow:K],ones(size(CVaR,1),1)*CVaR_star,'r-', 'LineWidth',1.5);
%     
%     t=sprintf('Case %d, b=%.2f, c=%.2f',experi,bEnlarge,cEnlarge);
%     title(t,'interpreter','latex','FontSize',15);
%     subplot(4,1,2);
%     plot(3*[1:skiprow:K],theta1,'b-', 'LineWidth',1.5);
%     ylabel('$\theta_k$','interpreter','latex','FontSize',15);
%     hold on
%     plot(3*[1:skiprow:K],theta2,'k-', 'LineWidth',1.5)
% %     plot([1:skiprow:K],theta3,'k-', 'LineWidth',1.5)
% %     plot([1:skiprow:K],theta4,'k-', 'LineWidth',1.5)
% %     plot([1:skiprow:K],theta5,'k-', 'LineWidth',1.5)
% %     plot([1:skiprow:K],theta6,'k-', 'LineWidth',1.5)
% %     plot([1:skiprow:K],theta7,'k-', 'LineWidth',1.5)
% %     plot([1:skiprow:K],theta8,'k-', 'LineWidth',1.5)
% %     plot([1:skiprow:K],theta9,'k-', 'LineWidth',1.5)
% %     plot([1:skiprow:K],theta10,'k-', 'LineWidth',1.5)
%     legend('$\theta^1_k$','$\theta^2_k$', 'interpreter','latex')
%     subplot(4,1,3);
%     plot([1:skiprow:K],dist_spco,'k-', 'LineWidth',1.5);
%     hold on
%     ylabel('$distance$','interpreter','latex','FontSize',15);
%     
%     
%     gcf
%     set(gcf, 'Position', [2 1 600 420]);
    
%     path_figure ='/Users/meichen/Documents/SPSA-Q/SPCO_all/SPCO/num_rate_sen_0527';
%     filename_figure_CVaR = sprintf('CVaR_sen%d_%.2f_%.2f',experi, bEnlarge, cEnlarge);
%     savefig(fullfile(path_figure, sprintf('%s.fig', filename_figure_CVaR)))
%     saveas(gcf,fullfile(path_figure, sprintf('%s.eps', filename_figure_CVaR)),'eps')
%     saveas(gcf,fullfile(path_figure, sprintf('%s.pdf', filename_figure_CVaR)),'pdf')
%     filename_CVaR_SPCO = sprintf('CVaR_SPCO_%d_%.2f_%d.csv', experi, alpha, d);
%     path_figuresave = '/Users/meichen/Documents/SPSA-Q/SPCO/numerical';
    % mkdir(path_figuresave);
%     saveas(gcf,fullfile(path_figuresave, sprintf('%s.eps', filename_CVaR_SPCO)),'eps')
%     saveas(gcf,fullfile(path_figuresave, sprintf('%s.pdf', filename_CVaR_SPCO)),'pdf')
%     
%     savefig(fullfile(path_figuresave, sprintf('%s.fig', filename_CVaR_SPCO)))
    
    