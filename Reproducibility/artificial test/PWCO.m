function [CVaR,index]=PWCO(K,alpha,experi,d,poolsize)%,skiprow) % 加distribution变量
    skiprow=1;
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
            case 1
            hk = h3(thetan,d);
            p = rand(nk_1,1);
            pn = cauchycdf(b_upper-h3(thetan,d)) *p + cauchycdf(b_lower-h3(thetan,d))*(1-p); %truncate
            X = tan(pi*(pn-0.5));
            Y=X+h3(thetan,d);
            Y_sort = sort(Y);
            D=h3_diff(thetan,d).*ones(1,nk_1);%d*nk_1
            q(k)=Y_sort(ceil((1-alpha)*nk_1));
            delta(:,k-1)=sum(D.*(Y>q(k)).',2)/(nk_1*alpha); %d*nk_1
            
            case 2
            p = rand(nk_1,1);
            X = -log(1-p); %nk *1
            hk = Ackley(thetan,d,alpha);
            Y = X * hk; %nk*1
            Y_sort = sort(Y);
            q(k)=Y_sort(ceil((1-alpha)*nk_1));
            D=Ackley_diff(thetan,d,alpha)*Y';
            delta(:,k-1)=sum(D.*(Y>q(k)).',2)/(nk_1*alpha); %d*nk_1
            case 3
            h5k=h5(thetan,d);
            h2k=h2(thetan,d);
            p = rand(nk_1,1);
            X = sqrt(2) * erfinv(2 * p - 1);
            Y = X * h2k + h5k; %nk*1
            Y_sort = sort(Y);
            q(k)=Y_sort(ceil((1-alpha)*nk_1));
            D=h2_diff(thetan,d)*Y'+h5_diff(thetan,d);
            delta(:,k-1)=sum(D.*(Y>q(k)).',2)/(nk_1*alpha);
            case 4
            p = rand(nk_1,1);
            X = sqrt(2) * erfinv(2 * p - 1);
            Y = Easom_x(thetan,d,X); %nk*1
            Y_sort = sort(Y);
            q(k)=Y_sort(ceil((1-alpha)*nk_1));
            D=Easom_x_diff(thetan,d,X);
            delta(:,k-1)=sum(D.*(Y>q(k)).',2)/(nk_1*alpha);
            
        end 
        thetatmp=theta(:,k-1)-ck*delta(:,k-1);
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
    plot(index,q,'k-.', 'LineWidth',1.5);
    % % plot([1:skiprow:K],q,'k-.', 'LineWidth',1.5);
    % hold on
          
      
    