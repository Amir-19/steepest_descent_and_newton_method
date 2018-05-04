% Question 7 part a
syms x1 x2 s;

%f1
f1 = x1^2+x2^2+4*x1-2*x2;

%for starting point (0,0)
iteration=steepest_descent(f1,x1,x2,s,[0,0],10^(-6));

result_string=sprintf('find the answer after %d number of iterations\n',iteration);
disp(result_string);

%for starting point (2,1.5)
iteration=steepest_descent(f1,x1,x2,s,[2,1.5],10^(-6));

result_string=sprintf('find the answer after %d number of iterations\n',iteration);
disp(result_string);

%f2
%for starting point (0,0)
f2 = 5*x1^2+x1^4-9*x1^2*x2+3*x2^2+2*x2^4+0.25*x1;

iteration=steepest_descent(f2,x1,x2,s,[0,0],10^(-6));

result_string=sprintf('find the answer after %d number of iterations\n',iteration);
disp(result_string);

%for starting point (2,1.5)
iteration=steepest_descent(f2,x1,x2,s,[2,1.5],10^(-6));

result_string=sprintf('find the answer after %d number of iterations\n',iteration);
disp(result_string);