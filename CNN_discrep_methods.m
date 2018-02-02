function [yp, resp, gamm, qgamm, itergamm]= CNN_discrep_methods(ys, B, sigma2, m,gamm)

if size(ys,2)>1
    ys=ys';
end

ns=length(ys);

% F matrix
if m==0
    F=eye(ns);
else
    c(1)=1;
    c(2,1)=-1;
    c(ns,1)=0;
    r(1)=1;
    r(1,ns)=0;
    % Delta random walk u(k)=u(k-1)+e
    DELTA=toeplitz(c,r);
    % F=delta^m
    F=DELTA;
    for k=2:m
        F=F*DELTA;
    end
end

% B from ar model
IB2=(ones(ns,1)./diag(B)).^(1/2); % B^(-1/2)
IB=diag(diag(B).^-1); %matrice
% init
conv=0; itergamm=0;

% ------------------------------------------------------------------------
% regualrization methods
gmin=1e-12;
gmax=1e+12;

if isempty(gamm)
    % Consistency 1
    while conv==0 & itergamm<=25
        itergamm=itergamm+1;
        gamm=10^((log10(gmin)+log10(gmax))/2);
        Dgamm=inv(IB+gamm*F'*F)*IB;
        uhat=Dgamm*ys;
        wess=sum((F*uhat).^2);
        qgamm=trace(Dgamm);
        lambda2=sigma2/gamm;
        if wess> lambda2*qgamm
            gmax= gamm; % lower gamm
        else
            gmin = gamm; % upper gamm
        end %if
        if abs((wess-lambda2*qgamm)/wess)<10e-06
            conv=1;
        end %if
    end %while
    
    % if not ....
    if conv==0
        itergamm=0;
        gmin=1e-10;
        gmax=1e+10;
    end
    % TRY Consistency 2
    while conv==0 & itergamm<=25
        itergamm=itergamm+1;
        gamm=10^((log10(gmin)+log10(gmax))/2);
        Dgamm=inv(IB+gamm*F'*F)*IB;
        uhat=Dgamm*ys;
        yp=uhat;
        resp=IB2.*(ys-yp);
        wrss=resp'*resp;
        qgamm=trace(Dgamm);
        if wrss> sigma2*(ns-qgamm)
            gmax= gamm; % lower gamm
        else
            gmin = gamm; % upper gamm
        end %if
        if abs((wrss-sigma2*(ns-qgamm))/wrss)<10e-06
            conv=1;
        end %if
    end
    
    % ELSE fixed lambda val
    if conv == 0
        gamm=5000;
        Dgamm=inv(IB+gamm*F'*F)*IB;
        uhat=Dgamm*ys;
        yp=uhat;
        resp=IB2.*(ys-yp);
        wrss=resp'*resp;
        qgamm=trace(Dgamm);
    else
        yp=uhat;
        resp=IB2.*(ys-yp);
        wrss=resp'*resp;
    end
else
    Dgamm=inv(IB+gamm*F'*F)*IB;
    uhat=Dgamm*ys;
    yp=uhat;
    resp=IB2.*(ys-yp);
    wrss=resp'*resp;
    qgamm=trace(Dgamm);
end


end