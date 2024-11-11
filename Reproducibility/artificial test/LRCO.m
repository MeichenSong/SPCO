function [CVaR,index]=LRCO(K,alpha,experi,d,poolsize)%,skiprow) % 加distribution变量

CVaR(1)=0;CVaR(2)=0;q(1)=0;q(2)=0;
% delta(:,1)=zeros(d,1);delta(:,2)=zeros(d,1);

% initialize theta
% switch funs
%     case 1
%          theta=4*rand(d,K)-2; c=1;
%     case 2
%         theta=([1:d]'+(4*rand(d,1)-2))*ones(1,K); c=0.5;theta_star = [1:d].';
%     case 3
% %         theta=4*rand(d,K)-2; c=10;
%         theta(:,1)=4*rand(d,1)-2;
%         theta(:,2)=4*rand(d,1)-2;B=d;theta_star = 0.5*[1:d].';
%     case 4
%         theta=3*rand(d,K)+1; c=0.75;
% end
switch experi
case 1 %h3
    theta = [1:d]'.*rand(d,2);
    theta_star = 0.5*[1:d].';
    theta_lower=zeros(d,1);
    theta_upper=[1:d].';

    b_upper = h3(zeros(d,1),d)+10000;
    b_lower = h3(theta_star,d)-10000;
    cdf_upper = cauchycdf(b_upper);
    cdf_lower = cauchycdf(b_lower);
    scale = cdf_upper - cdf_lower;
    transformed_level = scale*(1-alpha)+cdf_lower;
case 2
%         theta=zeros(d,K);
    theta = 4.*rand(d,2)-2;
    theta_lower = -2.*ones(d,1);
    theta_upper = 2.*ones(d,1);
    theta_star=zeros(d,1);
case 3
    theta = [1:d]'+rand(d,2)-0.5;
    theta_star = [1:d].';
    theta_upper = [1:d]'+0.5;
    theta_lower=[1:d]'-0.5;
case 4
        theta = pi.*rand(d,2)/2+3*pi/4;
%         theta_lower = (1*d-1)*pi/(1*d).*ones(d,1);
%         theta_upper = (1*d+1)*pi/(1*d).*ones(d,1);
%         theta_star = pi.*ones(d,1);
        theta_lower = -1.*ones(d,1);
        theta_upper = 1.*ones(d,1);
        theta_star = zeros(d,1);
end
% if alpha==0.05
%     b_upper = h3(-d*ones(d,1),d)+35;
%     b_lower = h3(theta_star,d)-35;
% elseif alpha==0.01
%     b_upper = h3(-d*ones(d,1),d)+35;
%     b_lower = h3(theta_star,d)-35;
% end

N=0;
index(1)=0;
index(2)=0;
for k=3:K
%     bk=bEnlarge*(1)/(k+10*bEnlarge);
%     ck=bEnlarge/(k+10*bEnlarge);
    ck=1/(k);
    nk_1=ceil((log(k-1))^4);
    N=N+nk_1;
    if N > poolsize
        break
    end
    index(k)=N;
    thetan=theta(:,k-1);
    switch experi
%          case 1
%              p = rand(nk_1,1);
%              X=sqrt(2) * erfinv(2 * p - 1);
%              Y=X+h3(theta(:,k-1),d);
%              Y_sort = sort(Y);
%              q(k)=Y_sort(ceil((1-alpha)*nk_1));
%              delta=((2*theta(:,k-1)-[1:d]')*((Y-h3(theta(:,k-1),d))')*diag(Y-q(k))*(Y>=q(k)))/(nk_1*alpha);
% %              p = rand(1,3);
% %              X=sqrt(2) * erfinv(2 * p - 1);
        case 1
        hk = h3(thetan,d);
        p = rand(nk_1,1);
%              X = tan(pi*(p-0.5));
%              Y=X+hk;
%              Y_sort = sort(Y);
%              q(k)=Y_sort(ceil((1-alpha)*nk_1));
%              delta=((2*thetan-[1:d]')*((2*(Y-hk)./(1+(hk-Y).^2))')*diag(Y-q(k))*(Y>=q(k)))/(nk_1*alpha);
        pn = cauchycdf(b_upper-h3(thetan,d)) *p + cauchycdf(b_lower-h3(thetan,d))*(1-p); %truncate
        X = tan(pi*(pn-0.5));
        Y=X+h3(thetan,d);
        Y_sort = sort(Y);
        q(k)=Y_sort(ceil((1-alpha)*nk_1));
        Fc = (1/((b_upper-hk)^2+1)-1/((b_lower-hk)^2+1))/(atan(b_upper-hk)-atan(b_lower-hk));
        delta=((2*thetan-[1:d]')*(Fc+(2*(Y-h3(thetan,d))./(1+(h3(thetan,d)-Y).^2))')*diag(Y-q(k))*(Y>=q(k)))/(nk_1*(alpha));
        case 2
        p = rand(nk_1,1);
        X = -log(1-p); %nk *1
        hk = Ackley(thetan,d,alpha);
        Y = X * hk; %nk*1
        Y_sort = sort(Y);
        q(k)=Y_sort(ceil((1-alpha)*nk_1));
        
        deri_Ackley = (4*exp(-0.2*thetan'*thetan/d)/((1-log(alpha))*d)).*thetan;
        deri_log_term1 = -deri_Ackley./hk;
        deri_log_term2 = deri_Ackley*Y'/(hk)^2; %d*nk
        deri_log = deri_log_term1+deri_log_term2;

        delta = (deri_log * diag(Y-q(k))*(Y>=q(k)))/(nk_1*(alpha));
        case 3
        h5k=h5(thetan,d);
        h2k=h2(thetan,d);
        p = rand(nk_1,1);
        X = sqrt(2) * erfinv(2 * p - 1);
        Y = X * h2k + h5k; %nk*1
        Y_sort = sort(Y);
        q(k)=Y_sort(ceil((1-alpha)*nk_1));
        
        deri_h5=prod(repmat(thetan-[1:d]',1,d).*(eye(d)==0)+eye(d),1); %1*d
        deri_h2=2*(thetan-[1:d]')/d;%d*1
        deri_log_term1 = -deri_h2/h2k;
        deri_log_term2 = -deri_h5'*(Y-h5k)'/h2k^2;
        deri_log_term3 = -deri_h2*((Y-h5k).^2)'/h2k^3;
        deri_log = deri_log_term1-deri_log_term2-deri_log_term3;
        delta = (deri_log * diag(Y-q(k))*(Y>=q(k)))/(nk_1*(alpha));
        case 4
        p = rand(nk_1,1);
        X = sqrt(2) * erfinv(2 * p - 1);
        Y = Easom_x(thetan,d,X); %nk*1
        Y_sort = sort(Y);
        q(k)=Y_sort(ceil((1-alpha)*nk_1));
        h6k = (thetan'*thetan)/d+0.1;
        
        deri_lognormal_miu = 2*thetan/(d*h6k);
        deri_log = deri_lognormal_miu*(log(Y)'-log(h6k)-1);
        % h6k = -prod(cos(thetan));
        % deri_h6 = prod(repmat(cos(thetan),1,d).*(eye(d)==0)+diag(sin(thetan)),1);%1*d
        % deri_lognormal_miu = deri_h6'/h6k-2*(thetan-pi);
        % deri_log = deri_lognormal_miu*(log(Y)'-logNormal_miu(thetan,d));
        delta = (deri_log * diag(Y-q(k))*(Y>=q(k)))/(nk_1*(alpha));
        
    end 
    thetatmp=theta(:,k-1)-ck*delta;
    thetatmp= thetatmp.*(thetatmp<=theta_upper)+theta_upper.*(thetatmp>theta_upper);
    theta(:,k)= thetatmp.*(thetatmp>=theta_lower)+theta_lower.*(thetatmp<theta_lower);
     
    hk=h3(theta(:,k),d);
    switch experi
    case 1
        centerk=h3(theta(:,k),d);
        scalek=1;
        true_q=cauchyinv(transformed_level,centerk);
        CVaR(k)= (scalek*log((b_upper-centerk)^2+scalek^2)+2*centerk*atan((b_upper-centerk)/scalek)/scalek-scalek*log((true_q-centerk)^2+scalek^2)-2*centerk*atan((true_q-centerk)/scalek)/scalek)/(2*pi*alpha*scale);
    case 2
        CVaR(k)=(-log(alpha)+1)*Ackley(theta(:,k),d,alpha);
    case 3
        hk_center = h5(theta(:,k),d);
        hk_scale = h2(theta(:,k),d);
        if alpha==0.05
            CVaR(k) = hk_scale * 2.062713+hk_center;
        elseif alpha==0.01
            CVaR(k) = hk_scale * 2.665214+hk_center;
        end
    case 4
        if alpha ==0.05
            CVaR(k) = 179.5081*0.5*Easom(theta(:,k),d);
        else
            CVaR(k) = 10262*0.5*Easom(theta(:,k),d);
        end
    end
end
    
switch experi
case 1
    switch d
    case 20
        switch alpha
            case 0.05
                t='$Case 1$, $\varphi = 0.95$, $d=20$';
            case 0.01
                t='$Case 1$, $\varphi = 0.99$, $d=20$';
        end
        
    case 50
        switch alpha
            case 0.05
                t='$Case 1$, $\varphi = 0.95$, $d=50$';
            case 0.01
                t='$Case 1$, $\varphi = 0.99$, $d=50$';
        end  
    end
    center_star = h3(theta_star,d);
    scale_star = 1;
    true_q=cauchyinv(transformed_level,center_star);
    CVaR_star= (scale_star*log((b_upper-center_star)^2+scale_star^2)+2*center_star*atan((b_upper-center_star)/scale_star)/scale_star-scale_star*log((true_q-center_star)^2+scale_star^2)-2*center_star*atan((true_q-center_star)/scale_star)/scale_star)/(2*pi*alpha*scale);
case 2
    switch d
    case 20
        switch alpha
            case 0.05
                t='Case 2, $\varphi = 0.95$, $d=20$';
            case 0.01
                t='Case 2, $\varphi = 0.99$, $d=20$';
        end
        
    case 50
        switch alpha
            case 0.05
                t='Case 2, $\varphi = 0.95$, $d=50$';
            case 0.01
                t='Case 2, $\varphi = 0.99$, $d=50$';
        end  
    end
    CVaR_star=(-log(alpha)+1)*Ackley(theta_star,d,alpha);
case 3
    switch d
    case 20
        switch alpha
            case 0.05
                t='Case 3, $\varphi = 0.95$, $d=20$';
            case 0.01
                t='Case 3, $\varphi = 0.99$, $d=20$';
        end
        
    case 50
        switch alpha
            case 0.05
                t='Case 3, $\varphi = 0.95$, $d=50$';
            case 0.01
                t='Case 3, $\varphi = 0.99$, $d=50$';
        end  
    end
    if alpha==0.05
        CVaR_star = h5(theta_star,d) + h2(theta_star,d) * 2.062713;
    elseif alpha==0.01
        CVaR_star = h5(theta_star,d) + h2(theta_star,d) * 2.665214;
    end
case 4
    switch d
    case 20
        switch alpha
            case 0.05
                t='Case 4, $\varphi = 0.95$, $d=2$';
            case 0.01
                t='Case 4, $\varphi = 0.99$, $d=2$';
        end
        
    case 50
        switch alpha
            case 0.05
                t='Case 4, $\varphi = 0.95$, $d=50$';
            case 0.01
                t='Case 4, $\varphi = 0.99$, $d=50$';
        end  
    end
    if alpha ==0.05
        CVaR_star = 179.5081*0.5*Easom(theta_star,d);
    else
        CVaR_star = 10262*0.5*Easom(theta_star,d);
    end
end

  
% CVaR = CVaR(1:skiprow:end);
% % q = q(1:skiprow:end);
% theta1 = theta(1,:);
% theta2 = theta(2,:);
% % theta1 = theta(1,1:skiprow:end);
% % theta2 = theta(2,1:skiprow:end);
% %     theta3 = theta(3,1:skiprow:end);
% %     theta4 = theta(4,1:skiprow:end);
% %     theta5 = theta(5,1:skiprow:end);
% %     theta6 = theta(6,1:skiprow:end);
% %     theta7 = theta(7,1:skiprow:end);
% %     theta8 = theta(8,1:skiprow:end);
% %     theta9 = theta(9,1:skiprow:end);
% % theta10 = theta(10,1:skiprow:end);
% subplot(3,1,1);
% % plot([1:skiprow:K],CVaR,'k-.', 'LineWidth',1.5);
% plot(index,CVaR,'k-.', 'LineWidth',1.5);
% ylabel('$CVaR(\theta_k)$','interpreter','latex','FontSize',15);
% hold on
% plot(index,ones(1,size(CVaR,2))*CVaR_star,'r-.', 'LineWidth',1.5);
% % plot([1:skiprow:K],ones(1,size(CVaR,2))*CVaR_star,'r-', 'LineWidth',1.5);
% title(t,'interpreter','latex','FontSize',15);
% subplot(3,1,2);
% plot(index,theta1,'b-.', 'LineWidth',1.5);
% % plot([1:skiprow:K],theta1,'b-.', 'LineWidth',1.5);
% ylabel('$\theta_k$','interpreter','latex','FontSize',15);
% hold on
% plot(index,theta2,'k-.', 'LineWidth',1.5)
% % plot([1:skiprow:K],theta2,'k-.', 'LineWidth',1.5)
%     % plot([1:skiprow:K],theta3,'k-', 'LineWidth',1.5)
%     % plot([1:skiprow:K],theta4,'k-', 'LineWidth',1.5)
%     % plot([1:skiprow:K],theta5,'k-', 'LineWidth',1.5)
%     % plot([1:skiprow:K],theta6,'k-', 'LineWidth',1.5)
%     % plot([1:skiprow:K],theta7,'k-', 'LineWidth',1.5)
%     % plot([1:skiprow:K],theta8,'k-', 'LineWidth',1.5)
%     % plot([1:skiprow:K],theta9,'k-', 'LineWidth',1.5)
%     % plot([1:skiprow:K],theta10,'k-', 'LineWidth',1.5)
% legend('$\theta^1_k$','$\theta^2_k$', 'interpreter','latex')
% subplot(3,1,3);
% plot(index,q,'k-.', 'LineWidth',1.5);
% % plot([1:skiprow:K],q,'k-.', 'LineWidth',1.5);
% hold on
      
  
