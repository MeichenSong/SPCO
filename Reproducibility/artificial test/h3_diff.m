function Y=h3_diff(theta,d) % d*1
    Y=2*theta-[1:1:d]';
    % Y=(theta-[1:1:d]')'*theta/(d^3);