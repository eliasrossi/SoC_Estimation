function c = costfcn2(eqdif,funch,x,xp,y,PAR,u)
% cost function

    if (PAR.mu == 0)
        c = [];
    else
        c = sqrt(PAR.mu).*(x-xp); % + norm(y(:,1)-channel(x,0))^2;
    end

    for t=1:PAR.N
        
        yp = funch(x,u(:,t),PAR);
        c =  [c ; y(:,t) - yp];

        x = eqdif(x,u(:,t),PAR);
    end

end