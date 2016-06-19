%%
 % Truncated conjugated gradient algorithm implementation
 % Steihaug and Toint algorithm
 % Ref.: Andrew R. Conn, Nicholas I. M. Gould, Philippe L. Toint,
 %   "Trust-Region Methods", SIAM, pp 202-207, 2000.
 % Ref : AMLET
 % Mai Anh Tien - maianhti@iro.uumontreal.ca
%%  
function step = btr_step_steihaug_toint(Op)    
    g =  Op.grad;
    v =  Op.grad;
    p = -Op.grad;
    Op.step = zeros(Op.n,1);
    normg0 = norm(g);
    k = 1;
    n = length(g);
    while norm(g) > normg0 * min(0.01, sqrt(normg0)) && k <= n
        Hp = Op.H * p;
        kappa = p'*Hp;
        if (kappa <= 0) 
            sigma = Op.step' * p;
            sqrnrm = p' * p;
            root = sqrt(sigma * sigma + sqrnrm * (Op.delta * Op.delta - Op.step' * Op.step));
            sigma = (root - sigma)/sqrnrm;
            Op.step = Op.step + sigma * p;
            k = n+1;
        else
            alpha = (g' * v)/kappa;
            sum = Op.step + alpha * p;
            if (sum' * sum >= Op.delta^2) 
                % leaving the trust region */
                sigma = Op.step' * p;
                sqrnrm = p' * p;
                root = sqrt(sigma * sigma + sqrnrm * (Op.delta^2 - Op.step' * Op.step));
                sigma = (root - sigma) / sqrnrm;
                Op.step = Op.step + sigma * p;
                k = n+1; % go out 
            else
                Op.step = sum;
                tmp = g' * v;
                % Compute g_k+1 
                g = g + alpha * Hp;
                v = g;
	            beta = (g' * v) / tmp;
	            p = beta * p - v;
            end            
        end;
        k = k+1;
    end
    step = Op.step;
end