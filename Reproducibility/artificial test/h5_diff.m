function Y=h5_diff(theta,d)
    for i = 1:d
        Y(i,1)=prod(theta-[1:d]')/(theta(i)-i);
    end
