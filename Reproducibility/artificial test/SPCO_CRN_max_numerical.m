function [CVaR]=SPCO_CRN_max_numerical(K,alpha,experi,d, WarmUp, bEnlarge, cEnlarge,skiprow,nk_2,empirical_indi) 

    CVaR=zeros(K,1);q=zeros(K,1);
    
    % initialize theta
    switch experi
    case 1 %h3
        theta = [1:d]'.*rand(d,K);
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
        theta = pi/2.*rand(d,K)-pi/4;
        theta_lower = -2.*ones(d,1);
        theta_upper = 2.*ones(d,1);
        theta_star=zeros(d,1);
    case 3
        theta = [1:d]'+rand(d,K)-0.5;
        theta_star = [1:d].';
        theta_upper = [1:d]'+0.5;
        theta_lower=[1:d]'-0.5;
    case 4
        theta = pi.*rand(d,K)/2+3*pi/4;
%         theta_lower = (1*d-1)*pi/(1*d).*ones(d,1);
%         theta_upper = (1*d+1)*pi/(1*d).*ones(d,1);
%         theta_star = pi.*ones(d,1);
        theta_lower = -1.*ones(d,1);
        theta_upper = 1.*ones(d,1);
        theta_star = zeros(d,1);
    end
    
    C=ceil(K*WarmUp);
    N=0;
    for k=2:K
        switch empirical_indi
            case 1
                nk_1=ceil((log(k))^4);
                N=N+nk_1;
            case 0
                nk_1=nk_2; N=N+nk_1;
        end
                
       
%         bk=bEnlarge*(2*C)/(k-1+C);
%         ck=cEnlarge*((2*C)^(0.167))/(k-1+C)^(0.167);
%         bk=bEnlarge*(1)/(k-1);
%         bk=bEnlarge*(1)/(k+10*bEnlarge);
%         ck=cEnlarge*(1)/(k+10*cEnlarge)^(0.167);
%         gammak=C/(k-1)^(5/9);
        
        bk=bEnlarge*(1)/((k+10*bEnlarge)^0.99);
        ck=cEnlarge*(1)/((k+10*cEnlarge)^0.1429);
        gammak=C/(k-1)^(0.5714);
        
        delta=2*(rand(d,1)>0.5)-1;
%         if experi == 4
%             delta=2*pi*(rand(d,1)>0.5)/d-pi/d;
%         end
%         
        thetan=theta(:,k-1);
        thetap=thetan+ck*delta;
        thetam=thetan-ck*delta;
        
        % sample input
        switch experi
        case 1
            p = rand(nk_1,1);
            
            centern=h3(thetan,d);centerp=h3(thetap,d);centerm=h3(thetam,d);
            scalen=1;scalep=1;scalem=1;
            l1n = (b_lower-centern)/scalen; l2n=(b_upper-centern)/scalen;
            l1p = (b_lower-centerp)/scalep; l2p=(b_upper-centerp)/scalep;
            l1m = (b_lower-centerm)/scalem; l2m=(b_upper-centerm)/scalem;
            pn = cauchycdf(l2n) *p(end,1) + cauchycdf(l1n)*(1-p(end,1));
            pn_q = cauchycdf(l2n) *p + cauchycdf(l1n)*(1-p);
            if size(p,1) == 1
                pp = cauchycdf(l2p) *p(1,1) + cauchycdf(l1p)*(1-p(1,1));
                pm = cauchycdf(l2m) *p(1,1) + cauchycdf(l1m)*(1-p(1,1));
            else
                pp = cauchycdf(l2p) *p(1:end-1,1) + cauchycdf(l1p)*(1-p(1:end-1,1));
                pm = cauchycdf(l2m) *p(1:end-1,1) + cauchycdf(l1m)*(1-p(1:end-1,1));
            end
            X_q_root=tan(pi*(pn-0.5));
            clear X
            X(:,1)=tan(pi*(pp-0.5));
            X(:,2)=tan(pi*(pm-0.5));
            
            X_q=tan(pi*(pn_q-0.5));
            
            
        case 2
            p = rand(nk_1,1);
            clear X
            if nk_1==1
                X=-log(1-p);
            else
                X = -log(1-p(1:end-1,1));
            end
            X_q_root = -log(1-p(end,1));
            X_q = -log(1-p);
        case 3
            p = rand(nk_1,1);
            clear X
            if nk_1==1
                X=sqrt(2) * erfinv(2 * p - 1);
            else
                X=sqrt(2) * erfinv(2 * p(1:end-1,1) - 1);
            end
            X_q_root = sqrt(2) * erfinv(2 * p(end,1) - 1);
            X_q=sqrt(2) * erfinv(2 * p - 1);
        case 4
            p = rand(nk_1,1);
            clear X
            if nk_1 ==1
                X=sqrt(2) * erfinv(2 * p - 1);
            else
                X=sqrt(2) * erfinv(2 * p(1:end-1,1) - 1);
            end
            X_q_root=sqrt(2) * erfinv(2 * p(end,1) - 1);
            X_q=sqrt(2) * erfinv(2 * p - 1);
        end

        switch experi
        case 1
            Y=X_q_root+h3(thetan,d); Yp=X(:,1)+h3(thetap,d); Ym=X(:,2)+h3(thetam,d); Y_q=X_q+h3(thetan,d);
        case 2
            Y=X_q_root*Ackley(thetan,d,alpha); Yp=X*Ackley(thetap,d,alpha); Ym=X*Ackley(thetam,d,alpha); Y_q=X_q*Ackley(thetan,d,alpha);
        case 3
            Y=X_q_root*h2(thetan,d)+h5(thetan,d); Yp=X*h2(thetap,d)+h5(thetap,d); Ym=X*h2(thetam,d)+h5(thetam,d); Y_q=X_q*h2(thetan,d)+h5(thetan,d);
        case 4
            Y=Easom_x(thetan,d,X_q_root); Yp=Easom_x(thetap,d,X); Ym=Easom_x(thetam,d,X); Y_q=Easom_x(thetan,d,X_q);
        end  
        
        switch empirical_indi
            case 1
                Y_sort = sort(Y_q);
                q(k)=Y_sort(ceil((1-alpha)*nk_1));
            case 0
                q(k)=q(k-1)+gammak*(1-alpha-(Y<=q(k-1)));
        end
        
        thetatmp = theta(:,k-1)-bk*(mean((max(q(k-1),Yp))-mean(max(q(k-1),Ym)))./(2*alpha*ck*delta));
%         q(k)=q(k-1)+gammak*(1-alpha-(Y<=q(k-1)));
%         indi = -(1-alpha-(Y<=q(k-1))<=0)+ (1-alpha-(Y<=q(k-1))>=0);
%         q(k)=q(k-1)+1*sign(1-alpha-(Y<=q(k-1)));
        
        thetatmp= thetatmp.*(thetatmp<=theta_upper)+theta_upper.*(thetatmp>theta_upper);
        theta(:,k)= thetatmp.*(thetatmp>=theta_lower)+theta_lower.*(thetatmp<theta_lower);
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
    
              
    CVaR = CVaR(1:skiprow:end);
    q = q(1:skiprow:end);
    theta1 = theta(1,1:skiprow:end);
    theta2 = theta(2,1:skiprow:end);
%     theta3 = theta(3,1:skiprow:end);
%     theta4 = theta(4,1:skiprow:end);
%     theta5 = theta(5,1:skiprow:end);
%     theta6 = theta(6,1:skiprow:end);
%     theta7 = theta(7,1:skiprow:end);
%     theta8 = theta(8,1:skiprow:end);
%     theta9 = theta(9,1:skiprow:end);
%     theta10 = theta(10,1:skiprow:end);
    subplot(3,1,1);
    plot([1:skiprow:K],CVaR,'k--', 'LineWidth',1.5);
    ylabel('$CVaR(\theta_k)$','interpreter','latex','FontSize',15);
    hold on
    plot([1:skiprow:K],ones(size(CVaR,1))*CVaR_star,'r-', 'LineWidth',1.5);
%     title(t,'interpreter','latex','FontSize',15);
    subplot(3,1,2);
    plot([1:skiprow:K],theta1,'b--', 'LineWidth',1.5);
    ylabel('$\theta_k$','interpreter','latex','FontSize',15);
    hold on
    plot([1:skiprow:K],theta2,'k--', 'LineWidth',1.5)
%     plot([1:skiprow:K],theta3,'k-', 'LineWidth',1.5)
%     plot([1:skiprow:K],theta4,'k-', 'LineWidth',1.5)
%     plot([1:skiprow:K],theta5,'k-', 'LineWidth',1.5)
%     plot([1:skiprow:K],theta6,'k-', 'LineWidth',1.5)
%     plot([1:skiprow:K],theta7,'k-', 'LineWidth',1.5)
%     plot([1:skiprow:K],theta8,'k-', 'LineWidth',1.5)
%     plot([1:skiprow:K],theta9,'k-', 'LineWidth',1.5)
%     plot([1:skiprow:K],theta10,'k-', 'LineWidth',1.5)
    legend('$\theta^1_k$','$\theta^2_k$', 'interpreter','latex')
    subplot(3,1,3);
    plot([1:skiprow:K],q,'k--', 'LineWidth',1.5);
    hold on
    ylabel('$q_k$','interpreter','latex','FontSize',15);
%     filename_CVaR_SPCO = sprintf('CVaR_SPCO_%d_%.2f_%d.csv', experi, alpha, d);
%     path_figuresave = '/Users/meichen/Documents/SPSA-Q/SPCO/numerical';
    % mkdir(path_figuresave);
%     saveas(gcf,fullfile(path_figuresave, sprintf('%s.eps', filename_CVaR_SPCO)),'eps')
%     saveas(gcf,fullfile(path_figuresave, sprintf('%s.pdf', filename_CVaR_SPCO)),'pdf')
%     
%     savefig(fullfile(path_figuresave, sprintf('%s.fig', filename_CVaR_SPCO)))
%     
    