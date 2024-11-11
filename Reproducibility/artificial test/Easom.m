function Y=Easom(theta,d)
    % term1=-prod(cos(theta).^2);
    % term2=exp(-(theta-pi)'*(theta-pi));
    % Y=(term1+1)*term2;

    term1=(theta'*theta)/d+0.1;
    term2=exp(1);
    Y=term1*term2;