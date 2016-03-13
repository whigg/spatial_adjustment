
tic

%format long;
format short;


X = [
    100.1030
    111.6010
    122.1810
    116.6920
    87.66100
    86.85500
    129.5510
    102.4480
    126.6760
    143.9770
    145.6870
    133.6100
    ];

Y = [
    100.0110
    109.0030
    144.0130
    168.0140
    134.1990
    106.2100
    161.8670
    90.16700
    96.81400
    115.7720
    140.4290
    163.0790
    ];

Distance = [
    14.5967
    49.2302
    69.9973
    36.5734
    59.2302
    24.6209
    ];

%Start Point

SP = [
    1
    1
    1
    2
    2
    3
    ];

%End Point

EP = [
    2
    3
    4
    3
    4
    4
    ];


%Directions Grad
DG =[
    0
    3.92813
    26.65768
    28.14169
    29.46106
    42.51533
    79.96743
    125.68087
    272.63544
    350.12485
    379.70302
    0
    28.95076
    114.46143
    170.86868
    205.16361
    233.14019
    236.90870
    239.06410
    252.26800
    306.12031
    365.52605
    0
    7.24369
    10.92332
    35.65047
    71.45038
    119.97324
    195.22762
    204.68251
    243.91905
    347.23966
    383.60808
    0
    3.71152
    9.75390
    24.10203
    29.54640
    45.87412
    66.81895
    86.84505
    97.16398
    370.06912
    388.29153
    0
    380.49429
    132.61652
    95.41328
    26.15275
    ];

%Directions Degrees DD:
DD = DG *(360/400); 

%Directions Radians DR:
%DR = degtorad(DD);
DR=[0.006004725
0.006004761
0.006005482
0.006002594
0.006004464
0.006008522
0.00600353
-0.060042837
0.006004595
0.0060013
0.006002864
-0.000885404
-0.000895906
-0.000886103
-0.000885849
-0.000893817
-0.000890511
-0.00088666
-0.000891942
-0.000899945
-0.000890824
0.008906963
-0.002637813
-0.002628415
-0.002627926
-0.002636812
-0.002643222
-0.002639072
-0.002638145
-0.002646034
-0.00263777
-0.002636193
0.026371401
-0.002419699
-0.002418774
-0.002418716
-0.002414115
-0.002410034
-0.002415644
-0.002415444
-0.002410641
-0.002413796
-0.002422643
0.024159505
-0.014280394
0.05714302
-0.014286201
-0.014289764
-0.014286661
];
%Stand Point
sTP = [
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
3
3
3
3
3
3
3
3
3
3
3
4
4
4
4
4
4
4
4
4
4
4
5
5
5
5
5
];

%Target Point
TP = [
    2
    11
    12
    3
    7
    4
    5
    6
    8
    9
    10
    1
    8
    9
    10
    11
    12
    7
    3
    4
    5
    6
    1
    8
    2
    9
    10
    11
    12
    7
    4
    5
    6
    1
    8
    2
    9
    3
    10
    11
    7
    12
    5
    6
    1
    6
    4
    3
    2
    ];

% if rx == zeros(12,1) & ry == zeros(12,1)
%     v = L - Ax*dx;
for loop=1:10
%while rx ~= zeros(12,1) & ry ~= zeros(12,1)
    X=X+rx;
    Y=Y+ry;
    % DISTANCE COMPUTATIONS

    %Approx distance ADs1:

    ADs1 = [];   %empty  ADs1

    for m=1:length(SP)
    i =SP(m);
    j =EP(m);
    
    dX = X(j) - X(i);            %compute and fill ADs1 matrix
    dY = Y(j) - Y(i);
    
    ADs1(m,1) = sqrt(dX^2 + dY^2);
    end

    %Distance L Matrix DsLM:

    DsLM = []; %empty DsLM

    for m=1:length(ADs1)
    i =Distance(m);
    j =ADs1(m);
    
    DsLM(m,1) = i-j;
    
    end

    %Distance A Matrix sA:

    sA = zeros(6,24); %empty sA

    for m=1:length(SP)
    i=SP(m);
    j=EP(m);
    s=ADs1(m);
    
    dX = X(j) - X(i);
    dY = Y(j) - Y(i);
    
    a=i*2;
    b=j*2;
    sA(m,a-1)   = -(dX/s);   %start point coefficients
    sA(m,a)     = -(dY/s);
    
    sA(m,b-1)   = (dX/s);   %end point coefficients
    sA(m,b)     = (dY/s);
    
    end

    %distance nuisance params sN:
    sN = zeros(6,5);
    for m=1:length(SP)
    i = SP(m);
    
    sN(m,i)=1;
    
    end


    %DIRECTIONS COMPUTATIONS

    %Approx Distances ADs2:

    ADs2 = [];   %empty  ADs2

    for m=1:length(sTP)
    i =sTP(m);
    j =TP(m);
    
    dX = X(j) - X(i);            %compute and fill ADs2 matrix
    dY = Y(j) - Y(i);
    
    ADs2(m,1) = sqrt(dX^2 + dY^2);
    end

    %Approx Directions ADr

    ADr = [];   %empty AD

    for m=1:length(sTP)
    i =sTP(m);
    j =TP(m);
    
    dX = X(j) - X(i);
    dY = Y(j) - Y(i);
    
    ADr(m,1) = atan(dY/dX);
    
    end

    %Direction L Matrix DrLM:

    DrLM = []; %empty DrLM

    for m=1:length(ADr)
    i = DR(m);
    j = ADr(m);
    
    DrLM(m,1) = i-j;
    
    end

    %Directions A Matrix dA:

    dA = zeros(49,24); %empty dA

    for m=1:length(sTP)
    i=sTP(m);
    j=TP(m);
    s=ADs2(m);
    
    dX = X(j) - X(i);
    dY = Y(j) - Y(i);
    
    a=i*2;
    b=j*2;
    dA(m,a-1)   = (dY/s);   %start point coefficients
    dA(m,a)     = -(dX/s);
    
    dA(m,b-1)   = -(dY/s);   %end point coefficients
    dA(m,b)     = (dX/s);
        end

    %direction nuisance params dN:
    dN = zeros(49,5);
    for m=1:length(sTP)
    i = sTP(m);
    
    dN(m,i)=1;
    
    end

    %CONCATENATIONS
    L = [DsLM;DrLM];

    % COMBINED A MATRIX 
    
    Ax = [sA;dA];  %observations A

    At = [sN;dN];

    %A = [sA,sN;dA,dN];

    A = [Ax At];

    %--------------------------------------------------------------------------

    %COMPUTATIONS:  PART II

    %WEIGHTS
    w1=1/(0.3^2);              %distance weights
    w2=1/(degtorad(1/3600)^2); %direction weights

    W = zeros(55);
    for m=1:6
       W(m,m)=w1;
    end
    for m=7:55
       W(m,m)=w2;
    end

    N = A'*W*A; %normal equation matrix of the form [Nxx, Nxt; Ntx, Ntx]

    Nxx = N(1:24,1:24);


    Nxt = N(1:24,25:end);
    Ntx = N(25:end,1:24);

    Ntt = N(25:end,25:end);

    %Reduced Normal Equation Matrix NXX: 2

    nxx = Nxx-Nxt*inv(Ntt)*Ntx ;

    %G MATRIX
    %Two dimensional Network with distances and directions as observations
    %%The scale datum defect is eliminated from distances observations.

    G = zeros(24,3);
    
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

    %Apply the G constraint matrix: 27*27

    NXX = [nxx,G;G',zeros(3,3)]; 

    QXX = inv(NXX);
    Qxx = QXX(1:24,1:24);

    %absolute vectors
    nx=Ax'*W*L;
    nt=At'*W*L;

    %reduced absolute vectors
    n_x = nx-Nxt*inv(Ntt)*nt;

    z=zeros(3,1);
    N_X = [nx;z];

    %PARAMETER MATRIX
    DX = Qxx*nx;
    %DX = Qxx*N_X;

    %params dx
    dx=DX(1:24);

    rx=[]; % X residuals
    ry=[]; % Y residuals

    for m=1:12    
    i=m*2-1;
    j=m*2;  
    
    rx(m,1)=dx(i); 
    ry(m,1)=dx(j);   
   
    end    
    
    %loop 
    %rx
    %ry
end
%%
%RESIDUAL MATRIX v
v = L - Ax*dx;

%Aposteriori Variance

apos = (v'*W*v)/(55-29);

%Covariance Matrix EXX
Exx = apos*Qxx;

%RESULTS
%%
%1.coordinates

coords =zeros(12,3);
for m=1:12
    coords(m,1)=m;
    coords(m,2)=X(m);
    coords(m,3)=Y(m);   
end
%coords

%%
%2. standard deviations
sig = diag(Exx);
sigXY = diag(Exx,1);
sigmas =zeros(12,3);

for m=1:12    
    i=m*2-1;
    j=m*2;
    sigX(m,1)=sig(i); 
    sigY(m,1)=sig(j);
    
    %display results
    sigmas(m,1)=m;
    sigmas(m,2)=sigX(m);
    sigmas(m,3)=sigY(m);  
   
end 

%sigmas

%%
%3. standard error ellipses
a=[];
b=[];
orient=[];
elements =zeros(12,4);
for m=1:12
    a(m)= 0.5*(sigX(m)+sigY(m))+sqrt(0.4*(sigX(m)-sigY(m))^2+sigXY(m));
    b(m)= 0.5*(sigX(m)+sigY(m))-sqrt(0.4*(sigX(m)-sigY(m))^2+sigXY(m));

    orient(m)= radtodeg( 0.5*atan((2*sigXY(m))/(sigX(m)-sigY(m))));

    %display results
    elements(m,1)=m;
    elements(m,2)=a(m);
    elements(m,3)=b(m);
    elements(m,4)=orient(m);
    
    %plot
    plot(X,Y)    

end
elements;

%%
% Global Model Test [ GMT ]
v1 = v(1:6);
v2 = v(7:55);

T1= (v1'*w1*v1)/(0.3^2);
T2= (v2'*w2*v2)/(degtorad(1/3600)^2);

%Chi-Square comparison at 95% C.I :
  %T = ;

toc
