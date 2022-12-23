function [] = heat_conduction(~)
%first non-linear equation.
err= input('saturation limit between two consecutive readings ');
gse= input('gauss-siedel saturation limit between two consecutive readings ');
c=input('max no. of iteration for gauss siedel method ');
l= input('length of the rod(in metres) '); 
n= input('number_of_points_of study '); %number of grid points
m= input('number_of_iterations '); %how many times you want to iterate, for achieving the steady state
k= input('thermal conductivity '); % thermal conductivity of a material
t= input('time of computation '); %for how long you want to compute
%c=k*t/m*(l/n);
p = input('lambda ');
%drichlet boundary conditions
flag= false;
flag2= false;
a=input('right hand side temp = ');
b=input('left hand side temp = ');
A = zeros(m+1,n+2); % the 2 D matrix was set to zero initial values, which was furthur used in Gauss-seidel method.
B = zeros(1,n);     %This is the comparator & substitution matrix.
C = zeros(1,n);     %This is the convergence check matrix.

% setting the boundary condition on the matrix element
for i=1:m+1
    A(i,1) = a;
end
for i=1:m+1
    A(i,n+2) = b;
end
disp(A);

%implicit calculation of the temprature profile, with an open loop method
i=1;
e =0;
o=0;
product=1;
while (i<= m )
    
j=2;
while( j<= n+1 ) 
       
     
    
    A(i+1,j)= (A(i,j) + c*(A(i+1,j+1)) +c*(A(i+1,j-1)) )/(1+2*c)  ; % this is the implicit equation.
    
    
    %if ((A(i+1,j))-(A(i,j)) < err ) convergence condition
     %   flag = true;
      
   j=j+1;
    
   end
    
   disp(A);
   disp('1st row of A'); % this is the first row of a after convergence
   
   
    j=2;   
    count=0;
    flag2= false;
    
    %
    while ( flag2== false)
       
        count=count+1;
       
      for k=1:n
        
      B(1,k)= (A(i,j) + c*(A(i+1,j+1)) +c*(A(i+1,j-1)) )/(1+2*c)  ;  % this is the implicit equation
    
      j=j+1;    
      end
      
      j=2;
    disp(B);
    disp('new A as B');
    disp(C);
    disp('c');
    
    for s=2:n+1
        e = A(i+1,s) - B(1,s-1);
        %disp(e);
        %disp('e');
        %disp('iteration number');
        if (gse < e )
            C(1,s-1)= 1;
        end
        
    end
    
    disp(C);
    disp('C matrix');
        
        for o=1:n
            product=C(1,o)*product;
        end
        
        %here we see that the saturation has been reached.
    if (product == 1 || count== c)
        flag2 = true;
        disp('gs saturation reached for this row');
        
    end
    
    
        for s=2:n+1
            A(i+1,s)=B(1,s-1); % substitution matrix
        end
        
    end
   
i=i+1;
disp(i);   
end


%if (flag == true)
%        disp( 'steady state reached ' );    
    


A(m+1,n+2)= b;

disp(A);
disp('final');

hold on
for i=1:m+1
plot(A(i,:))

end

end
