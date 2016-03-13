
tic

%format long;
format short;


X = load('X_Coordinates.csv');     %Coordinates
Y = load('Y_Coordinates.csv');   
D = load('Distances.csv');         %Distances

SP = load('Start_Points.csv');     %Start Point

EP = load('End_Points.csv');       %End Point

DG =load('Directions_Grads.csv');  %Directions Grad

DD = DG *(360/400);                %Directions Degrees DD:

DR = degtorad(DD);                 %Directions Radians DR:

sTP = load('Stand_Point.csv');     %Stand Point

TP = load('Target_Point.csv');     %Target Point



%COMPUTATIONS:

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

while W ~= ones(55,55)             % OUTLIER DETECTION LOOP
               
    %%
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
    DsLM = []; 

     for m=1:length(ADs1)
        i =D(m);
        j =ADs1(m);
    
        DsLM(m,1) = i-j;
    
     end

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
     ADs2 = [];   

     %Approx Directions ADr
     ADr = [];   

     for m=1:length(sTP)
        i =sTP(m);
        j =TP(m);
        
        dX = X(j) - X(i);           
        dY = Y(j) - Y(i);
    
        ADs2(m,1) = sqrt(dX^2 + dY^2);

        ADr(m,1) = atan(dY/dX);
     end
    

     %Direction L Matrix DrLM:
     DrLM = []; 

     for m=1:length(ADr)
        i = DR(m);
        j = ADr(m);
    
        DrLM(m,1) = i-j;
    
     end

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
        dA(m,a-1)   = (dY/s);                 %start point coefficients
        dA(m,a)     = -(dX/s);
    
        dA(m,b-1)   = -(dY/s);                %end point coefficients
        dA(m,b)     = (dX/s);

        dN(m,i)=1;                            %direction nuisance params dN:
    
     end


     %CONCATENATIONS
    
     %combined L matrix:
     L = [DsLM;DrLM];

     %combined A matrix:     
     Ax = [sA;dA];                            %observations A

     At = [sN;dN];                            %nuisances A

     
     A = [Ax At];
    
     %%

     N = A'*W*A;                              %normal equation matrix of the form [Nxx, Nxt; Ntx, Ntx]

     Nxx = N(1:24,1:24);

     Nxt = N(1:24,25:end);
     Ntx = N(25:end,1:24);                    %extracted components

     Ntt = N(25:end,25:end);     

     nxx = Nxx-Nxt*inv(Ntt)*Ntx ;

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

     QXX = inv(NXX);                          %Cofactor matrix of parameters Qxx: (24*24)
     Qxx = QXX(1:24,1:24);

     %absolute vectors
     nx=Ax'*W*L;
     nt=At'*W*L;

     %reduced absolute vectors
     n_x = nx-Nxt*inv(Ntt)*nt;

     %PARAMETER MATRIX

     dx=Qxx*nx;                               %params dx

     rx=[];                                   % X params
     ry=[];                                   % Y params

     for m=1:12    
        i=m*2-1;
        j=m*2;  
    
        rx(m,1)=dx(i); 
        ry(m,1)=dx(j);  
   
     end 

     %RESIDUAL MATRIX 

     v = L - Ax*dx ;    
     
     outlier = max(abs(v))                    %maximum of residuals
      
     [r,c] = find(v==max(v(:)))              %index of max residual
                          
     W(r,r) = 0;
     
     
    
end

%%

% apos = (v'*W*v)/(55-24);                      %Aposteriori Variance

% Exx = apos*Qxx;                               %Covariance Matrix EXX

% %RESULTS

% %%
% coords =zeros(12,3);                          %1.coordinates

% for m=1:12
%     coords(m,1)=m;
%     coords(m,2)=X(m);
%     coords(m,3)=Y(m);   
% end

% %%
% sig = diag(Exx);                              %2. standard deviations
% sigXY = diag(Exx,1);
% sigmas =zeros(12,3);

% for m=1:12    
%     i=m*2-1;
%     j=m*2;
%     sigX(m,1)=sig(i); 
%     sigY(m,1)=sig(j);    
    
%     %display results
%     sigmas(m,1)=m;            
%     sigmas(m,2)=sigX(m);
%     sigmas(m,3)=sigY(m);  
   
% end 

% %%


% a=[];
% b=[];
% orient=[];
% elements =zeros(12,4);                        %3. standard error ellipses


% for m=1:12
%     a(m)= (0.5*(sigX(m)+sigY(m))+sqrt(0.4*(sigX(m)-sigY(m))^2+sigXY(m)))/1000;
%     b(m)= (0.5*(sigX(m)+sigY(m))-sqrt(0.4*(sigX(m)-sigY(m))^2+sigXY(m)))/1000;

%     orient(m)= radtodeg( 0.5*atan((2*sigXY(m))/(sigX(m)-sigY(m))));

%     %display results
%     elements(m,1)=m;
%     elements(m,2)=abs(a(m));
%     elements(m,3)=abs(b(m));
%     elements(m,4)=orient(m);
   
% end 

%  x= [];
%  y= [];

% %plot network
% for m=1:length(sTP)    
%     a=sTP(m);
%     b=TP(m);

%     x(m,1) = X(a);
%     x(m,2) = X(b);
%     y(m,1) = Y(a);
%     y(m,2) = Y(b);
    
%     plot(y,x,'-r');

%     hold on    

% end

%  p= [];
%  q= [];
% for m=1:length(elements)

%     t=-2*pi:0.01:2*pi;

%     p=X(m)+(elements(m,2)*sin(t)); 
%     q=Y(m)+(elements(m,3)*cos(t));  

%     ellipse = plot(q,p,'-g'); 

%     % dr=[0 0 1];
%     % origin=[p q];
%     % rotate(plot(q,p,'-g'),dr,elements(m,4),origin);   

%     text(Y(m),X(m),num2str(m));

%     hold on  
            
% end
% xlabel('Y Coordinates');
% ylabel('X Coordinates');
% title('Adjusted Network');



% %%

% % GLOBAL MODEL TEST [ GMT ]:

% T = (v'*W*v)/1;    %test hypothesis

% %Chi-Square comparison at 95% C.I :
% chi = 46.979;

% if chi <= T
%     display('Accept the null hypothesis')
%     else
%         display('Reject the null hypothesis')
% end

toc

