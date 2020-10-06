

function [z, Sol]=CostFunction(x)

    n=numel(x);

    f1=x(1);
    
    g=1+9/(n-1)*sum(x(2:end));
    
    h=1-sqrt(f1/g);
    
    f2=g*h;
    
    f3=1/(f2+f1);
    
    z=[f1
       f2
       f3];
   Sol.f1=f1;
   Sol.f2=f2;
   Sol.f3=f3;
   Sol.x=x;
end