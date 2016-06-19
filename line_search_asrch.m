% More and Thuente Arc Search
%
% This is a modified implementation of the algorithm documented in:
%
% J. J. More and D. J. Thuente. Line Search Algorithms with Guaranteed
% Sufficient Decrease. TOMS 20-3. September 1994. pg 286-307.
%
% It attempts to find a step length stp such that:
%
% Sufficient decrease is met:
% f(x+s) <= f(x) + ftol*(ginit*stp + 0.5*min(ncur,0)*stp^2)
%
% Curvature condition is met:
% |g(x+s)| <= gtol*|ginit + min(ncur,0)*stp|
%
% s = arc(stp) is the displacement vector.
% ginit the first derivative along the search arc at stp=0.
%
% It makes the assumption that ftol <= gtol and thus does not require the
% modified interval updating rules. It also uses bisection to compute the
% next trial step instead of polynomial interpolation.
%
% This version uses a function of a vector and is efficient with
% evaluations.
%
% Input:
% fcn = a handle to a function that returns the value and first
% derivative. Usage in code: [f g] = fcn(x)
% x = current point
% f = function value at x
% g = gradient at x
% arc = search arc function. Usage in code: [s ds] = arc(stp)
% stp = initial step length
% ncur = negative curvature parameter
% ftol = sufficient decrease parameter, mu in the paper
% gtol = curvature condition parameter, eta in the paper
% xtol = the algorithm terminates if the width of the interval is less
% than xtol.
% stpmin = minimum step length allowed. It is acceptable to set stpmin=0.
% stpmax = maximum step length allowed.
% maxfev = maximum number of function evaluations
% fid = file identifier for test output (optional)
%
% Output:
% x = final point
% f = final function value
% g = final gradient
% stp = step length
% info = termination flag
% nfev = number of function evaluations.
%
% Termination Flags:
% info = 1 if stp satisfies the descent and curvature condition
% info = 2 if interval size is less than xtol
% info = 3 if algorithm has exceeded maxfev
% info = 4 if stpmin > 0 and stp == stpmin
% info = 5 if stp == stpmax
% info = 6 if stp == stpmax & strong wolfe conditions met
% info = 7 if rounding errors prevent progress
%%
function [x f g stp info nfev] = line_search_asrch(fcn,x,f,g,arc,stp,...
                                         ncur,ftol,gtol,xtol,...
                                         stpmin,stpmax,maxfev,...
                                         fid,bisect)

  % list of variables and parameters
  % extrap = parameter for extrapolations
  % bracket = true once bracket containing solution is found
  % info = status flag for output
  % nfev = number of function evaluations
  % s = displacement vector from arc
  % ds = derivative of displacement vector from arc
  %
  % stx = step size at "l" point
  % fx = function value at "l" point
  % dx = derivative of search function at "l" point
  %
  % sty = step size at "u" point
  % fy = function value at "u" point
  % dy = derivative of search function at "u" point
  %
  % stp = trial step size
  % fp = function value at trial step size
  % dp = derivative of search function at trial step size
  %
  % mfx = modified function value at "l" point
  % mdx = modified derivative value at "l" point
  %
  % mfp = modified function value at trial point
  % mdp = modified derivative value at trial point
  %
  % Note al and au define the bounds of the bracket if one is found. al and au
  % are the endpoints of the bracket, but are not ordered.
  %
  % finit = initial function value
  % ginit = initial gradient
  % amin = minimium step size
  % amax = maximum step size
  % ucase = arc search update case

  % optional input
  if nargin < 15 || isempty(bisect)
    bisect = 0;
  end
  
  % parameters
  xtrapu = 4;
  p66 = 0.66;

  % flags
  bracket = false;
  info = 0;

  % counters
  nfev = 0;

  % interval width tracker
  width = stpmax - stpmin;
  width1 = 2*width;
  
  % inital values
  [s ds] = arc(0);

  stx = 0;
  fx = f;
  dx = g'*ds;

  finit = fx;
  ginit = dx;

  fp = 0;
  dp = 0;

  sty = 0;
  fy = 0;
  dy = 0;

  % formatting & printing
  if nargin < 14 || isempty(fid) || fid <= 0
    print_flag = false;
  else
    print_flag = true;
  end

  ha_str = '%4s %1s %12s %12s %12s %12s %12s|%4s\n';
  d1_str = '%4d %1d %12g %12g %12g %12g %12g';
  d2_str = '|%4s\n';

  if print_flag
    fprintf(fid,ha_str,'nfev','b','stx','sty','stp','fp','dp','case');
  end

  while(1)

    % set the bounds of the interval
    if bracket
      stmin = min(stx,sty);
      stmax = max(stx,sty);
    else
      stmin = stx;
      %stmin = stp + xtrapl*(stp-stx);
      stmax = stp + xtrapu*(stp-stx);
    end

    % safeguard the trial step size
    stp = max(stp,stpmin);
    stp = min(stp,stpmax);
    
    % If an unusual termination is to occur then let
    % stp be the lowest point obtained so far.
    if ((bracket && (stp <= stmin || stp >= stmax)) ...
        || nfev >= maxfev-1 ...
        || (bracket && stmax-stmin <= xtol*stmax))
      stp = stx;
    end
    
    % evaluate function
    [s ds] = arc(stp);
    [f g] = fcn(x+s);
    fp = f;
    dp = g'*ds;
    nfev = nfev + 1;

    if print_flag
      fprintf(fid,d1_str,nfev,bracket,stx,sty,stp,fp,dp);
    end
    
    % compute modified function values
    mstx = stx;
    mfx = fx - finit - ftol*(ginit*stx + 0.5*min(ncur,0)*stx^2);
    mdx = dx - ftol*(ginit + min(ncur,0)*stx);
    
    %mstp = stp;
    mfp = fp - finit - ftol*(ginit*stp + 0.5*min(ncur,0)*stp^2);
    mdp = dp - ftol*(ginit + min(ncur,0)*stp);
    
    msty = sty;
    mfy = fy - finit - ftol*(ginit*sty + 0.5*min(ncur,0)*sty^2);
    mdy = dy - ftol*(ginit + min(ncur,0)*sty);
    
    % convergence tests
    
    % terminate if rounding errors prevent progress
    if bracket && (stp <= stmin || stp >= stmax)
      info = 7;
    end
    
    % terminate at stpmax
    if stp == stpmax && mfp <= 0 && mdp < 0
      info = 5;
    end
    
    % terminate at stpmin
    if stpmin > 0 && stp == stpmin && (mfp > 0 || mdp >= 0)
      info = 4;
    end
    
    % terminate if interval is too small
    if bracket && (stmax - stmin < xtol*stmax)
      info = 2;
    end
    
    % terminate if reached maximum number of function evaluations
    if nfev >= maxfev
      info = 3;
    end
    
    % terminate if strong wolfe conditions are met
    if fp <= finit + ftol*(ginit*stp + 0.5*min(ncur,0)*stp^2) ...
             && abs(dp) <= gtol*abs(ginit + min(ncur,0)*stp)
      info = 1;
    end
    
    % if strong wolfe conditions are met with at == stpmax
    if info == 1 && stp == stpmax
      info = 6;
    end
    
    if info
      x = x + s;
      
      if print_flag
        fprintf(fid,d2_str,['t' num2str(info)]);
      end
      
      return
    end
    
    % update the interval
    if mfp > mfx
      % case U1
      %stx = stx; fx = fx; dx = dx;
      sty = stp; fy = fp; dy = dp;
      bracket = true;
      ucase = 1;
    elseif mfp <= mfx && mdp*(stx-stp) > 0
      % case U2
      stx = stp; fx = fp; dx = dp;
      %sty = sty; fy = fy; dy = dy;
      ucase = 2;
    else % mfp <= mfx && mdp*(stx-stp) < 0
      % case U3
      sty = stx; fy = fx; dy = dx;
      stx = stp; fx = fp; dx = dp;
      bracket = true;
      ucase = 3;
    end
    
    % print the case
    if print_flag
      fprintf(fid,d2_str,['u' num2str(ucase)]);
    end
    
    % compute new trial step size
    if bracket && bisect
      % bisect if desired
      stp = stx + 0.5*(sty - stx);
    else
      % compute new step using interpolation
      stp = line_search_astep(mstx,mfx,mdx,msty,mfy,mdy,stp,mfp,mdp,bracket,stmin,stmax);
    end
    
    % safeguard the step and update the interval width tracker
    if bracket
      if (abs(sty-stx) >= p66*width1)
        stp = stx + 0.5*(sty - stx);
      end
      width1 = width;
      width = abs(sty-stx);
    end
    
  end

end