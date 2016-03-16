tic

%format long;
format short;


X = load('X_Coordinates.csv');     %Coordinates
Y = load('Y_Coordinates.csv');   
D = load('Error_distances.csv');         %Distances

SP = load('Start_Points.csv');     %Start Point

EP = load('End_Points.csv');       %End Point

DG =load('Directions_Grads.csv');  %Directions Grad

DD = DG *(360/400);                %Directions Degrees DD:

dr = degtorad(DD);                 %Directions Radians DR:

sTP = load('Stand_Point.csv');     %Stand Point

TP = load('Target_Point.csv');     %Target Point

dxTest = load('dxTest.csv');        %LoopFunction :: parameter matrix


  

%COMPUTATIONS:  PART II

%WEIGHTS
w1=1/(0.0003^2);              %distance weights
w2=1/(degtorad(1/3600)^2); %direction weights

W = zeros(55);
for m=1:6
    W(m,m)=w1;
end
for m=7:55
    W(m,m)=w2;
end


%%
%while dxTest ~= zeros(24,1);
for m=1:5
    
    % DISTANCE COMPUTATIONS

    %Approx distance ADs1:
    ADs1 = [];   

    for m=1:length(SP)
       i =SP(m);
       j =EP(m);

       dX = X(j) - X(i);            
       dY = Y(j) - Y(i);

       ADs1(m,1) = sqrt(dX^2 + dY^2);
    end

    %Distance L Matrix DsLM:
    DsLM = D - ADs1;    

     %Distance A Matrix sA:
     sA = zeros(6,24);           

     %distance nuisance params sN:
     sN = zeros(6,5);

     for m=1:length(SP)
        i=SP(m);
        j=EP(m);
        s=ADs1(m);
    
        dX = X(j) - X(i);
        dY = Y(j) - Y(i);
    
        a=i*2;
        b=j*2;
        sA(m,a-1)   = -(dX/s);                %start point coefficients
        sA(m,a)     = -(dY/s);
    
        sA(m,b-1)   = (dX/s);                 %end point coefficients
        sA(m,b)     = (dY/s);
         
        sN(m,i)=1;                            %distance nuisance params sN:
    
     end


     %%
     %DIRECTIONS COMPUTATIONS

     %Approx Distances ADs2:
     ADs2 = zeros(49,1);   

     %Approx Directions ADr
     ADr = zeros(49,1);   

     for m=1:length(sTP)
        i =sTP(m);
        j =TP(m);
        
        dX = X(j) - X(i);           
        dY = Y(j) - Y(i);
    
        ADs2(m,1) = sqrt(dX^2 + dY^2);

        if ((dX>0) &&(dY>0));
          ADr(m,1) = atan (abs(dY/dX));
         elseif ((dX<0) && (dY>0));
          ADr(m,1) = pi -(atan (abs(dY/dX)));
         elseif ((dX<0) && (dY<0));
          ADr(m,1) = atan (abs(dY/dX)) + pi;
         else
          ADr(m,1) = (2*pi)-atan (abs(dY/dX));    
        end
        
         %ADr(m,1) = atan(dY/dX);
     end
     
     %create and apply a swing to the directions observations
       for m=1:11;j=1;
        swing(m,j) = atan((Y(2,1)-Y(1,1))/(X(2,1)-X(1,1)));
        m=12:22;j=1;
        swing(m,j) = atan((Y(1,1)-Y(2,1))/(X(1,1)-X(2,1))) + pi;
        m = 23:33;j=1;
        swing(m,j) = atan((Y(3,1)-Y(1,1))/(X(3,1)-X(1,1))) + pi;
        m= 34:44;
        swing(m,j) = atan((Y(4,1)-Y(1,1))/(X(4,1)-X(1,1))) + pi;
        m = 45:49 ;
        swing(m,j) = atan((Y(5,1)-Y(1,1))/(X(5,1)-X(1,1))) ; 
        Dr= dr+swing; 
       end 
       
      for m=1:length(Dr) 
        if (Dr(m)>0 & Dr(m)<=pi);
          Obs_DR(m) = Dr(m);

        elseif (Dr(m)<=1.5*pi & Dr(m)>pi );
         Obs_DR(m) =  Dr(m);
     
        elseif (Dr(m)<=2*pi  & Dr(m)>1.5*pi);
           Obs_DR(m) = Dr(m);
    
        elseif (Dr(m)>2*pi);
          Obs_DR(m) = Dr(m)-2*pi;
    
        else (Dr(m)<=0);
         Obs_DR(m) = Dr(m) + 2*pi;
     
        end
       end
  

     %Direction L Matrix DrLM:
        
     DR = Obs_DR';
     
     DrLM = DR-ADr;
     
     %Directions A Matrix dA:
     dA = zeros(49,24);

     %direction nuisance params dN:
     dN = zeros(49,5); 

     for m=1:length(sTP)
        i=sTP(m);
        j=TP(m);
        s=ADs2(m);
    
        dX = X(j) - X(i);
        dY = Y(j) - Y(i);
    
        a=i*2;
        b=j*2;
        dA(m,a-1)   = (dY/s^2);                 %start point coefficients
        dA(m,a)     = -(dX/s^2);
    
        dA(m,b-1)   = -(dY/s^2);                %end point coefficients
        dA(m,b)     = (dX/s^2);

        dN(m,i)=1;                            %direction nuisance params dN:
    
     end


     %CONCATENATIONS
    
     %combined L matrix:
     L = [DsLM;DrLM];

     %combined A matrix:     
     Ax = [sA;dA];                            %observations A

     %At = [sN;dN];                            %nuisances A
      At = [zeros(6,5);dN];
     
     A = [Ax At];
    
     %%

     N = A'*W*A;                              %normal equation matrix of the form [Nxx, Nxt; Ntx, Ntx]
     
     Nxx = N(1:24,1:24);

     Nxt = N(1:24,25:end);
     Ntx = N(25:end,1:24);                    %extracted components

     Ntt = N(25:end,25:end);
     
     
     %absolute vectors
     nx=Ax'*W*L;
     nt=At'*W*L;

     %reduced absolute vectors
     n_x = nx-(Nxt*(inv(Ntt))*nt);
     
     nxx = Nxx-(Nxt*(inv(Ntt))*Ntx);

     %G matrix

     G = zeros(24,3);

     %Two dimensional Network with distances and directions as observations
     %The scale datum defect is eliminated from distances observations.
    
     for m=1:12
    
        n=m*2;
        a=1; b=2; c=3;
    
        G(n-1,a) = 1;
        G(n-1,b) = 0;
        G(n-1,c) = Y(m);
    
        G(n,a) = 0;
        G(n,b) = 1;
        G(n,c) = -X(m);
    
     end

     %Apply the G constraint matrix: (27*27)

     NXX = [nxx,G;G',zeros(3,3)];             %Reduced Normal Equation Matrix NXX: (27*27) 

     QXX = inv(NXX);                       %Cofactor matrix of parameters Qxx: (24*24)
         
     n_xx=[n_x ;zeros(3,1)];

     %PARAMETER MATRIX
     
     DX = QXX*n_xx ;                   %params dx
      
     rx=[];                                   % X params
     ry=[];                                   % Y params

     dx = DX(1:24);
     for m=1:12    
        i=m*2-1;
        j=m*2;  
    
        rx(m,1)= dx(i); %DX(i);
        ry(m,1)= dx(j); %DX(j);
    end  
     

    %RESIDUAL MATRIX 

     v = L - Ax*dx ;    
     
     % outlier = max(abs(v)) ;                   %maximum of residuals
      
     % [r,c] = find(v==max(v(:)));              %index of max residual
                          
     % W(r,r) = 0;                              %change weights at outlier index
     
    meanResidual = sum(v)/55;
    sqrMean= power((v-meanResidual),2);               %compute mean standard deviation for the data
    standard_dev = sqrt(sum(sqrMean)/55);

    p = zeros(length(v),1);

    for m= 1:length(v)
      if (v(m) <= 0.5*standard_dev)
        p(m) = 1;                              %determine multiplying factors to vary the weight
     
      else
        p(m) = exp(-v(m));
      end
      %set = 0.5*standard_dev
      %check = v(m)
      W(m)*p(m);
    end
     
end

%%

%RESIDUAL MATRIX 
%dx1=[DX;zeros(2,1)];

%v = L - (A*dx1); 
%v = L - Ax*dx; 

%GMM=(v'*W*v)
apos = (v'*W*v)/(55-24);                      %Aposteriori Variance

Qxx = QXX(1:24,1:24);
Exx = apos*Qxx;                               %Covariance Matrix EXX

%RESULTS

%%
coords =zeros(12,2);                          %1.coordinates

for m=1:12
    %coords(m,1)=m;
    coords(m,1)=X(m);
    coords(m,2)=Y(m);   
end

%%
var = diag(Exx);                              %2. standard deviations
covXY = diag(Exx,1);
sigmas =zeros(12,2);

for m=1:12    
    i=m*2-1;
    j=m*2;
    sigX(m,1)=sqrt(var(i)); 
    sigY(m,1)=sqrt(var(j));    
    
    %display results
    %sigmas(m,1)=m;            
    sigmas(m,1)=sigX(m);
    sigmas(m,2)=sigY(m);  
    sigmas(m,3)=covXY(m*2-1);
   
end 

%%

a=[];
b=[];
orient= zeros(12,1);
elements =zeros(12,3);                        %3. standard error ellipses


for m=1:12
    a(m)= sqrt((0.5*(sigX(m)^2+sigY(m)^2)+sqrt(0.4*(sigX(m)^2-sigY(m)^2)^2+covXY(m))));
    b(m)= sqrt((0.5*(sigX(m)^2+sigY(m)^2)-sqrt(0.4*(sigX(m)^2-sigY(m)^2)^2+covXY(m))));


    orient(m)= radtodeg( 0.5*atan((2*covXY(m))/(sigX(m)^2-sigY(m)^2)));

    %display results
    %elements(m,1)=m;
    elements(m,1)=abs(a(m));
    elements(m,2)=abs(b(m));
    elements(m,3)=orient(m);   
end 

 x= [];
 y= [];

%plot grid
for i = 1:12;
    x=X(i);
    y=Y(i);
    grid on;
    hold on;
    plot(x,y);
end
xlabel('X Coordinates');
ylabel('Y Coordinates');
title('Adjusted Network');

%plot network

 for m=1:length(sTP)    
     a=sTP(m);
     b=TP(m);

      p = X(a);
      g = X(b);
      q = Y(a);
      r = Y(b);

      x = [p g];
      y = [q r];
      
     plot(y,x,'-r');
     hold on   

 end

 p= [];
 q= [];
for m=1:length(elements)
 
    t=-2*pi:0.01:2*pi;

    p=X(m)+(elements(m,1)*sin(t)*25); 
    q=Y(m)+(elements(m,2)*cos(t)*25);  

    ellipse = plot(q,p,'-g'); 

%     % dr=[0 0 1];
%     % origin=[p q];
%     % rotate(plot(q,p,'-g'),dr,elements(m,4),origin);   

    text(Y(m),X(m),num2str(m));

    hold on  
            
end




%%

% GLOBAL MODEL TEST [ GMT ]:

T = (v'*W*v)         %test hypothesis

%Chi-Square comparison at 95% C.I :
chiU = 47.6;                    %set a check against lower and upper chi limits
chiL = 16.168;

if (T>chiL && T < chiU)
    display('Accept the null hypothesis')
else
    display('Reject the null hypothesis')
end

toc



