function Y=Easom_x_diff(theta,d,x)
    % x: n*1
    % theta: d*1
    % term1=-prod(cos(theta).^2);
    % term2=exp(x-(theta-pi)'*(theta-pi));
    % Y=(term1+1)*term2;
    term1=2*theta/d;
    term2=exp(x+1);
    Y=term1*term2';% d*n