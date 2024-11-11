function Y=Easom_x(theta,d,x)
    % term1=-prod(cos(theta).^2);
    % term2=exp(x-(theta-pi)'*(theta-pi));
    % Y=(term1+1)*term2;
    term1=(theta'*theta)/d+0.1;
    term2=exp(x+1);
    Y=term1*term2;