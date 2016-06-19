
% This function computes a safeguarded step for a search
% procedure and updates an interval that contains a step that
% satisfies a sufficient decrease and a curvature condition.
%
% The parameter stx contains the step with the least function
% value. If brackt is set to true (1) then a minimizer has
% been bracketed in an interval with endpoints stx and sty.
% The parameter stp contains the current step.
% The subroutine assumes that if brackt is set to true then
%
% min(stx,sty) < stp < max(stx,sty),
%
% and that the derivative at stx is negative in the direction
% of the step.
%
% The subroutine statement is
%
% stf = line_search_astep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,stpmin,stpmax)
%
% where
%
% stx is a double precision variable.
% On entry stx is the best step obtained so far and is an
% endpoint of the interval that contains the minimizer.
%
% fx is a double precision variable.
% On entry fx is the function at stx.
%
% dx is a double precision variable.
% On entry dx is the derivative of the function at
% stx. The derivative must be negative in the direction of
% the step, that is, dx and stp - stx must have opposite
% signs.
%
% sty is a double precision variable.
% On entry sty is the second endpoint of the interval that
% contains the minimizer.
%
% fy is a double precision variable.
% On entry fy is the function at sty.
%
% dy is a double precision variable.
% On entry dy is the derivative of the function at sty.
%
% stp is a double precision variable.
% On entry stp is the current step. If brackt is set to true
% then on input stp must be between stx and sty.
%
% fp is a double precision variable.
% On entry fp is the function at stp
%
% dp is a double precision variable.
% On entry dp is the the derivative of the function at stp.
%
% brackt is an logical variable.
% On entry brackt specifies if a minimizer has been bracketed.
% Initially brackt must be set to .false.
%
% stpmin is a double precision variable.
% On entry stpmin is a lower bound for the step.
%
% stpmax is a double precision variable.
% On entry stpmax is an upper bound for the step.
%
% MINPACK-1 Project. June 1983
% Argonne National Laboratory.
% Jorge J. More' and David J. Thuente.
%
% MINPACK-2 Project. November 1993.
% Argonne National Laboratory and University of Minnesota.
% Brett M. Averick and Jorge J. More'.
%%
function stpf = line_search_astep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,stpmin,stpmax)

% parameter
  p66 = 0.66;
  
  sgnd = dp*(dx/abs(dx));

  if (fp > fx)
    % First case: A higher function value. The minimum is bracketed.
    % If the cubic step is closer to stx than the quadratic step, the
    % cubic step is taken, otherwise the average of the cubic and
    % quadratic steps is taken.

    theta = 3.0*(fx-fp)/(stp-stx) + dx + dp;
    s = max([abs(theta) abs(dx) abs(dp)]);
    gamma = s*sqrt((theta/s)^2-(dx/s)*(dp/s));
    if (stp < stx)
      gamma = -gamma;
    end
    p = (gamma-dx) + theta;
    q = ((gamma-dx)+gamma) + dp;
    r = p/q;
    stpc = stx + r*(stp-stx);
    stpq = stx + ((dx/((fx-fp)/(stp-stx)+dx))/2.0)*(stp-stx);
    if (abs(stpc-stx) < abs(stpq-stx))
      stpf = stpc;
    else
      stpf = stpc + (stpq-stpc)/2.0;
    end
    %brackt = true;

  elseif (sgnd < 0.0)
    % Second case: A lower function value and derivatives of opposite
    % sign. The minimum is bracketed. If the cubic step is farther from
    % stp than the secant step, the cubic step is taken, otherwise the
    % secant step is taken.

    theta = 3.0*(fx-fp)/(stp-stx) + dx + dp;
    s = max([abs(theta) abs(dx) abs(dp)]);
    gamma = s*sqrt((theta/s)^2-(dx/s)*(dp/s));
    if (stp > stx)
      gamma = -gamma;
    end
    p = (gamma-dp) + theta;
    q = ((gamma-dp)+gamma) + dx;
    r = p/q;
    stpc = stp + r*(stx-stp);
    stpq = stp + (dp/(dp-dx))*(stx-stp);
    if (abs(stpc-stp) > abs(stpq-stp))
      stpf = stpc;
    else
      stpf = stpq;
    end
    %brackt = true;

  elseif (abs(dp) < abs(dx))
    % Third case: A lower function value, derivatives of the same sign,
    % and the magnitude of the derivative decreases.

    % The cubic step is computed only if the cubic tends to infinity
    % in the direction of the step or if the minimum of the cubic
    % is beyond stp. Otherwise the cubic step is defined to be the
    % secant step.

    theta = 3.0*(fx-fp)/(stp-stx) + dx + dp;
    s = max([abs(theta) abs(dx) abs(dp)]);

    % The case gamma = 0 only arises if the cubic does not tend
    % to infinity in the direction of the step.

    gamma = s*sqrt(max(0.0,(theta/s)^2-(dx/s)*(dp/s)));
    if (stp > stx)
      gamma = -gamma;
    end
    p = (gamma-dp) + theta;
    q = (gamma+(dx-dp)) + gamma;
    r = p/q;
    if (r < 0.0 && gamma ~= 0.0)
      stpc = stp + r*(stx-stp);
    elseif (stp > stx)
      stpc = stpmax;
    else
      stpc = stpmin;
    end
    stpq = stp + (dp/(dp-dx))*(stx-stp);

    if (brackt)

      % A minimizer has been bracketed. If the cubic step is
      % closer to stp than the secant step, the cubic step is
      % taken, otherwise the secant step is taken.

      if (abs(stpc-stp) < abs(stpq-stp))
        stpf = stpc;
      else
        stpf = stpq;
      end
      if (stp > stx)
        stpf = min(stp+p66*(sty-stp),stpf);
      else
        stpf = max(stp+p66*(sty-stp),stpf);
      end
    else

      % A minimizer has not been bracketed. If the cubic step is
      % farther from stp than the secant step, the cubic step is
      % taken, otherwise the secant step is taken.

      if (abs(stpc-stp) > abs(stpq-stp))
        stpf = stpc;
      else
        stpf = stpq;
      end
      stpf = min(stpmax,stpf);
      stpf = max(stpmin,stpf);
    end

  else
    % Fourth case: A lower function value, derivatives of the same sign,
    % and the magnitude of the derivative does not decrease. If the
    % minimum is not bracketed, the step is either stpmin or stpmax,
    % otherwise the cubic step is taken.

    if (brackt)
      theta = 3.0*(fp-fy)/(sty-stp) + dy + dp;
      s = max([abs(theta) abs(dy) abs(dp)]);
      gamma = s*sqrt((theta/s)^2-(dy/s)*(dp/s));
      if (stp > sty)
        gamma = -gamma;
      end
      p = (gamma-dp) + theta;
      q = ((gamma-dp)+gamma) + dy;
      r = p/q;
      stpc = stp + r*(sty-stp);
      stpf = stpc;
    elseif (stp > stx)
      stpf = stpmax;
    else
      stpf = stpmin;
    end
  end

end