%   Accept candidate ...
function isAccepted = btr_accept_candidate(rho)
    global Op;
    isAccepted = (rho >= Op.ETA1);
end