
function [] = btr_update_radius(rho)
    global Op;
    normstep = norm(Op.step);
    if (rho >= Op.ETA2)
        Op.delta = min(4 * normstep , 2* Op.delta);
    else
        if rho <= Op.ETA1
            Op.delta = Op.delta * 0.5;
        else
            Op.delta = Op.delta * 0.7;
        end
    end    
end