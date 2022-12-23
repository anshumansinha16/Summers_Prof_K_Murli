function [] = coupled_problem(~)



%first non-linear equation.

err= input('saturation limit between two consecutive readings ');
gse= input('gauss-siedel saturation limit between two consecutive readings ');
c=input('max no. of iteration for gauss siedel method ');

l= input('length of the rod(in metres) '); 
n= input('number_of_points_of study '); %number of grid points
m= input('number_of_iterations '); %how many times you want to iterate, for achieving the steady state
k= input('thermal conductivity '); % thermal conductivity of a material
ti= input('time of computation '); %for how long you want to compute

a=input('thermal diffusitivity ');

pn=input('pn ');
q=input('q ');
ep=input('ep ');
ko=input('ko ');
lu=input('lu ');
biq=input('biq ');
bim=input('bim ');
bims= bim*(1-(1-ep)*pn*ko*lu);
Ts=input('Surrounding temp ');
To=input('Initial temp ');
Q= (q*l)/(k*(Ts-To));
%c=k*t/m*(l/n);

%drichlet boundary conditions
flag= false;
flag2= false;
%t=(a*ti)/((l*l));
t=input('t ');

x=input('inpput x ');

% T & M are used for the initial zero matrix , which will be later used for
% gs method.
T = zeros(m,n);
M = zeros(m,n);
A = zeros(1,n);
B = zeros(1,n);
C = zeros(1,n);

i=1;
e1 =0;
e2 =0;
o=0;
product=1;

while (i<= m-1 )

    
j=2;

% initial conditions

T(i,1)= (4*T(i,2) - T(i,3) + 2*x*Q )/3 ;
M(i,1)= (M(i,3) -4*M(i,2) -pn*Q*2*x)/3;


% disretised initial conditions

T(i,n)=(biq+ 4*(T(i,n-1) - T(i,n-2))/(2*x) - (1-ep)*ko*lu*bim*(1 - M(i,n)))/ (3/(2*x) + biq);
M(i,n)=(bims+ (-4*(M(i,n-1) + T(i,n-2)))/2*x - pn*biq + - T(i,n)*pn*biq )/ (3/(2*x) + bims);
    

   while( j<= n-1 && i<=m-1 )
      
     T(i+1,1)= (4*T(i+1,2) - T(i+1,3) + 2*x*Q )/3 ;  
     M(i+1,1)= (M(i+1,3) -4*M(i+1,2) -pn*Q*2*x)/3;
         
    
     % discretised equations
    
    T(i+1,j)= (((T(i+1,j+1)+ T(i+1,j-1))/(x*x) ) + (T(i,j)/t) -(ep*ko*(M(i+1,j))) - (((M(i,j))/t))) / ((1/t) + 1/(x*x) )  ;
    M(i+1,j)= (((lu*(M(i+1,j+1)+ M(i+1,j-1))/(x*x) )) - ((lu*pn)*( T(i+1,j-1) - T(i+1,j) + T(i+1,j+1)) / (x*x)) + ((M(i,j))/ t ) ) / (1/t + lu/x*x ) ;
    
    
    %discretised boudary conditions
     
     T(i+1,n)=(biq+ 4*(T(i,n-1) - T(i,n-2))/(2*x) - ((1-ep)*ko*lu*bim*(1 - M(i,n))))/(3/(2*x) + biq);
     
     M(i+1,n)=(bims+ (-4*(M(i,n-1) + T(i,n-2)))/2*x - pn*biq + - T(i,n)*pn*biq )/(3/(2*x) + bims);
     
     
     %if ((A(i+1,j))-(A(i,j)) < err )
     %   flag = true;
      
      
    
    
    j=j+1;
    
   end
   
  
   disp(T);
   disp(M);
   disp('1st row of T & M');
   
   
    j=2;
    count=0;
    flag2= false;
    while ( flag2== false)
       
        count=count+1;
       
      
      A(1,1)= (M(i+1,3) -4*M(i+1,2) -pn*Q*2*x)/3;
      B(1,1)= (4*T(i+1,2) - T(i+1,3) + 2*x*Q )/3 ;
      
      for k=2:n-1
          
         
      %if( j~= n)
      A(1,k)= (((T(i+1,j+1)+ T(i+1,j-1))/(x*x) ) + (T(i,j)/t) -(ep*ko*(M(i+1,j))) - (((M(i,j))/t))) / ((1/t) + 1/(x*x) ) ;
               
      B(1,k)= (((lu*(M(i+1,j+1)+ M(i+1,j-1))/(x*x) )) - ((lu*pn)*( T(i+1,j-1) - T(i+1,j) + T(i+1,j+1)) / (x*x)) + ((M(i,j))/ t ) ) / (1/t + lu/x*x ) ;
      %end

      
      j=j+1;
    
      end
      
            
      A(1,n)=(biq+ 4*(T(i+1,n-1) - T(i+1,n-2))/(2*x) - (1-ep)*ko*lu*bim*(1- M(i+1,n)))/ (3/(2*x) + biq);
      B(1,n)=(bims+ (-4*(M(i+1,n-1) + T(i+1,n-2)))/2*x - pn*biq + - T(i+1,n)*pn*biq )/ (3/(2*x) + bims);
      
      
      
      j=2;
 
    disp(A);
    disp(B);
    disp('new A and B of T and M');
    disp(C);
    disp('c');
    
    for s=1:n
        e1 = T(i+1,s) - A(1,s);
        e2 = M(i+1,s) - B(1,s);
        %disp(e);
        %disp('e');
        %disp('iteration number');
        if (gse < e1 && gse < e2 )
            C(1,s)= 1;
        end
     end
    
    disp(C);
    disp('C matrix');
        
        for o=1:n
            product=C(1,o)*product;
        end
        
    if (product == 1 || count== c)
        flag2 = true;
        disp('gs saturation reached for this row'); % convergence of the row elemnts
        
    end
    
%if (flag == true)
%        disp( 'steady state reached ' );       
        for s=1:n
            T(i+1,s)=A(1,s); % substitution matrix
            M(i+1,s)=B(1,s);
        end
        
        disp(T);
        disp(M);
        disp('new T & M after replacement with A & B');
        
    end
   
i=i+1; % changing rows, changing time step
disp(i);   
end
disp(T);
disp(M);
disp('final');
plot(T(:,n))

end
