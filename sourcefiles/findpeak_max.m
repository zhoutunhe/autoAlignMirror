function [xpeak, ypeak, max_f] = findpeak_max(f,subpixel)
%FINDPEAK Find extremum of matrix.
%   [XPEAK,YPEAK,MAX_F] = FINDPEAK(F,SUBPIXEL) finds the extremum of F,
%   MAX_F, and its location (XPEAK, YPEAK). F is a matrix. MAX_F is the maximum
%   absolute value of F, or an estimate of the extremum if a subpixel
%   extremum is requested.
%
%   SUBPIXEL is a boolean that controls if FINDPEAK attempts to estimate the
%   extremum location to subpixel precision. If SUBPIXEL is false, FINDPEAK
%   returns the coordinates of the maximum absolute value of F and MAX_F is
%   max(abs(F(:))). If SUBPIXEL is true, FINDPEAK fits a 2nd order
%   polynomial to the 9 points surrounding the maximum absolute value of
%   F. In this case, MAX_F is the absolute value of the polynomial evaluated
%   at its extremum.
%
%   Note: Even if SUBPIXEL is true, there are some cases that result
%   in FINDPEAK returning the coordinates of the maximum absolute value
%   of F:
%   * When the maximum absolute value of F is on the edge of matrix F.
%   * When the coordinates of the estimated polynomial extremum would fall
%     outside the coordinates of the points used to constrain the estimate.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision $  $Date: 2004/10/20 17:54:47 $

% get absolute peak pixel

 [max_f, imax] = max(f(:));  %changed to actual value
[ypeak, xpeak] = ind2sub(size(f),imax(1));
    
if ~subpixel
    return % return absolute peak
elseif subpixel == 1    %use polynomial fitting
     if xpeak==1 %move away from edge
        xpeak = 2;
    elseif xpeak==size(f,2)
        xpeak = size(f,2)-1;
    end
    if ypeak==1
        ypeak = 2;
    elseif ypeak==size(f,1)
        ypeak = size(f,1)-1;
    end   
    
    % fit a 2nd order polynomial to 9 points  
    % using 9 pixels centered on irow,jcol    
    u = f(ypeak-1:ypeak+1, xpeak-1:xpeak+1);
    u = u(:);
    x = [-1 -1 -1  0  0  0  1  1  1]';
    y = [-1  0  1 -1  0  1 -1  0  1]';    

    % u(x,y) = A(1) + A(2)*x + A(3)*y + A(4)*x*y + A(5)*x^2 + A(6)*y^2
    X = [ones(9,1),  x,  y,  x.*y,  x.^2,  y.^2];
    
    % u = X*A
    A = X\u;

    % get absolute maximum, where du/dx = du/dy = 0
    x_offset = (-A(3)*A(4)+2*A(6)*A(2)) / (A(4)^2-4*A(5)*A(6));
    y_offset = -1 / ( A(4)^2-4*A(5)*A(6))*(A(4)*A(2)-2*A(5)*A(3));
    
  
    
    xpeak = xpeak + x_offset;
    ypeak = ypeak + y_offset;    
    
    % Calculate extremum of fitted function
    max_f = [1 x_offset y_offset x_offset*y_offset x_offset^2 y_offset^2] * A;
    
elseif subpixel == 2   % use interpolation

    if xpeak<3 %move away from edge
        xpeak = 3;
    elseif xpeak>size(f,2)-2
        xpeak = size(f,2)-2;
    end
    if ypeak<3
        ypeak = 3;
    elseif ypeak>size(f,1)-2
        ypeak = size(f,1)-2;
    end
    


    % use interpolation to find the minimum value  
    % using 25 pixels centered on irow,jcol    
    u = f(ypeak-2:ypeak+2, xpeak-2:xpeak+2);
    u_interp = interp2(-2:2,-2:2,u,-2:0.001:2,[-2:0.001:2]','cubic');
    [max_f, umin] = min(u_interp(:));
    [y_offset, x_offset] = ind2sub(size(u_interp),umin(1));
    y_offset = (y_offset-2001)/1000;
    x_offset = (x_offset-2001)/1000;

    
    xpeak = xpeak + x_offset;
    ypeak = ypeak + y_offset;    
    
    % Calculate extremum of fitted function

    
end
