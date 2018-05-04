function   iteration = newton(f,x1,x2,s,start_point,thereshold)
%
%   finds the minimum of multivariate function 'f(x1,x2)' with newton
%   approach from starting point of 'start_point' and the
%   termination condition having the accurecy of thereshold for striaght
%   line
%
%   f  -             symbolic function
%   x1 -             first variate
%   x2 -             second variate
%   s  -             parameter s of steepest descent
%   start_point -    example [0 0]
%   thereshold  -    termination accurecy
%
%   
    %keep track of how many iteration will it take to calculate
    iteration = 0;
    
    %calculate the gradiant of f
    grad_f = [diff(f,x1) diff(f,x2)];
    grad_2f = [diff(grad_f(1),x1)   diff(grad_f(1),x2)
             diff(grad_f(2),x1)   diff(grad_f(2),x2)];
    
    %calculate the delta for finding the step and accurecy
    delta = subs(grad_f,[x1,x2],[start_point(1),start_point(2)]);
    delta_2 = subs(grad_2f,[x1,x2],start_point);
    
    %assign the value of starting point to the current point
    current_point = start_point;
    
    %only for the first iteration
    step = 1;
    while norm(delta)*step > thereshold
        
        iteration = iteration + 1;
        
        %calculate the next point point based on delta, s and current point
        x_next = [current_point(1),current_point(2)] - s .* (inv(delta_2) * delta')';
        
        %symbolic substitution to have the value at next point
        fval = subs(f,[x1,x2],[x_next(1),x_next(2)]);
        
        %calculate the step size
        step = double( solve(  diff(fval,s) ,s)  );
        step = abs(step(1));
        
        %update current point
        current_point = double([current_point(1),current_point(2)] - step * (inv(delta_2) * delta')');
        
        %calculate the next step delta by symbolic substitution
        delta = subs(grad_f,[x1,x2],[current_point(1),current_point(2)]);
        
        %calculate objective value with symbolic substitution
        f_value = double(subs(f,[x1,x2],[current_point(1),current_point(2)]));
        
        
        %printing the progress
        result_string=sprintf('iteration=%d, x1=%.6f, x2=%.6f, step=%.6f f(x1,x2)=%.6f',...
        iteration,current_point(1),current_point(2),step,f_value);
        disp(result_string);
        
    end

end