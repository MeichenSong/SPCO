function Y=h1_rate(theta,d)
    Y=((theta-[1:1:d]')'*(theta-[1:1:d]')+1)/d;