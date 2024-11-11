function Y=Ackley(theta,d,alpha) %d*1
    % term1=-20*exp(-0.2*sqrt(theta'*theta/d));
    % term2=-exp(mean(cos(2*pi*theta)));
    % Y=(term1 + term2 +21+exp(1))/(1-log(alpha));
    term1=-10*exp(-0.2*theta'*theta/d)+11;
    % term1=-1*exp(-0.2*theta'*theta/d)+2;
    Y=term1/(1-log(alpha));