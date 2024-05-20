%-------------------------------------------------------------------------%
% Test Optimal sampling - with Pre-processing
% Description: We solve the weighted-LSP approximation of a function on an
% irregular domain, from samples points under a optimal measure.
% Programer: Juan Manuel Cardenas
% Date: July 16 - 2019 / Last modification: July 16 - 2019
%-------------------------------------------------------------------------%
 
%% Set up
tic
clear all; close all; clc;
space=' ';
parpool;                        % Start parallelize 
index_type = 1;                 % Index set HC
d_values = [2 3 5 10];%[15];%   % Dimension
m = 0;                          % number of k_points: M = m*N or M=N*log(N)
domain = -1;                    % 0 cube\corner, 1 cube\sphere, -1 annular  
f = 1;                          % type of function: f=1 , f= 0 Genz
Trials = 10;                    % number of trials
rad = 1/2;                        % annular's radio
N_max= 1000;                    % biggest size of |HC|
max_div = 10;                   % nmax with d=15 is 6, with d = 10 and 10>d is 20.
 
%% Run
 
Error_box = [];                     %zeros(Trials,length(d_values),max_div);
Min_sv = [];                        %zeros(Trials,length(d_values),max_div);
MIN_SV = [];
ERROR_BOX = [];
M_VALUES = [];
N_VALUES = [];
n_values = zeros(length(d_values),max_div);
N_values = zeros(length(d_values),max_div);
M_values = zeros(length(d_values),max_div);
 
for p = 1:length(d_values)
      
    d = d_values(p);
    %% Generate K points i.i.d from unit measure 
    K = 20*N_max;
    Z = [];
    if domain == 0                      % cube\ corner
        while length(Z) < K
         
            x = 2*rand(1,d)-rand(1,d);
         
            if norm(x,inf) <= 1
                if sum(x)<= 1  
                    Z = [Z;x];
                    disp(num2str(length(Z)));
                end
            end
         
        end
    elseif domain == 1                  % cube\ sphere rad = 1/2
        while length(Z) < K
         
             [15];% x = 2*rand(1,d)-rand(1,d);
         
            if norm(x,inf) <= 1
                if  norm(x) >= 1/2 
                    Z = [Z;x];
                    disp(num2str(length(Z)));
                end
            end
         
        end
    else                                  % annular             
        while length(Z) <= K
         
            u = randn(1,d);
            x = rad*(u/norm(u))*rand^(1/d);
         
            if norm(x)>rad/4
                Z = [Z ; x];
            end
         
        end
    end
    Z = Z(1:K,:);
     
    %% Pre - processing
 
    %--- Find nmax and Nmax ---%
     
    HCsize = 0;
    n_0 = 0;
     
    while HCsize < N_max
        n_0 = n_0 + 1;
        HCsize = length(HC_index(n_0,d));
    end
     
    nmax = n_0 -1;
    Nmax = length(HC_index(nmax,d));
    n_values(p,:) = round(linspace(1,nmax,max_div));
     
    %--- Re-order the index set ---%
     
    for i = 1 : length(n_values)
         
        I = HC_index(n_values(p,i),d);
        N_values(p,i) = length(I);
         
        if i == 1
            J = I;
        else
            C = setdiff(I',J','rows');              % Index in I\J
            C = C';                                   
            J = [J  C];                             % adding K to J indices
        end
    end
 
   %--- Generate B_max matrix ---% 
                                      
    B = zeros(K,Nmax);                              % legendre pol. evaluate in grid 
 
    for i = 1:K                                     % Loop over the points
        z = Z(i,:);                     
        L = LegMat(z',nmax+1);                      % evaluate legendre pol. in point z
        for j = 1:Nmax                              % Loop over the indexs
            Lij = zeros(d,1); 
            for k = 1:d                             % Loop over the components of indexs
                Lij(k,1) = L(k,J(k,j)+1); 
            end
            B(i,j) = prod(Lij);                     % Tensor product-type
        end
    end
     
    B_max = sqrt(1/K)*B;                            % Montecarlo factor            
    [Q,R] = qr(B_max,0);                            % Reduce QR-factorization
 
    %--- Measure mu ---%
     
    MU = abs(Q.^2);                                
  
    %--- Function right-hand side ---%
     
    if f == 1
       f_grid = func(Z);
    else
       f_grid = F_Genz(Z);
    end
     
%%  Compute W-LSP approximation 
 
    %--- M points ---%
     
    if m == 0                          
        M_values = round(N_values(p,:).*log(N_values(p,:)));  
    else                
        M_values = m*N_values(p,:);
    end
      
    for l = 1:length(N_values(p,:))
         
        N = N_values(p,l);
        M = M_values(l);
        mu = (1/N)*sum(MU(:,1:N)');
         
        parfor t = 1:Trials
             
            %--- Choosing M points i.i.d under d\mu  ---%
         
            I = datasample(1:K,M,'Replace',true,'Weights',mu);
 
            %--- Data matrix ---%
             
            W = 1./mu(I);
            D = diag(sqrt(W));                   % Save weights associate to M pts 
            A = (sqrt(K)/sqrt(M))*D*Q(I,1:N);        
             
            %--- Right-hand side ---%
         
            if f == 1 
                F = f_grid(I);
            end
            b = (1/sqrt(M))*D*F;
 
            %--- Compute least-squares approximation ---%
         
            c = A\b;
            s_min(t) = min(svd(A));                    % save min singular value
         
            %--- Compute error ---%
         
            f_tilde = sqrt(K)*Q(:,1:N)*c;
            Error(t) = norm(f_grid - f_tilde);     
        end
         
        %--- save ---%
         
        Min_sv(:,l) = s_min;
        Error_box(:,l) = Error;
         
        disp(' ');
        disp(['--------------------------------- DONE ---------------------------------']);
        disp(['dimension = ',num2str(d),space,'|| N = ',num2str(N),space,'M = ',num2str(M),space,'|| Iteration = ',num2str(l),space,'|| Measure = First ']);
        disp(['Err(1) = ',num2str(Error(1)),space,'|| s_min(1) = ',num2str(s_min(1))]);
        disp(['Err(2) = ',num2str(Error(2)),space,'|| s_min(2) = ',num2str(s_min(2))]);
        disp(['Err(3) = ',num2str(Error(3)),space,'|| s_min(3) = ',num2str(s_min(3))]);
        disp(['Err(4) = ',num2str(Error(4)),space,'|| s_min(4) = ',num2str(s_min(4))]);
        disp(['Err(5) = ',num2str(Error(5)),space,'|| s_min(5) = ',num2str(s_min(5))]);
        disp(['Err(6) = ',num2str(Error(6)),space,'|| s_min(6) = ',num2str(s_min(6))]);
        disp(['Err(7) = ',num2str(Error(7)),space,'|| s_min(7) = ',num2str(s_min(7))]);
        disp(['Err(8) = ',num2str(Error(8)),space,'|| s_min(8) = ',num2str(s_min(8))]);
        disp(['Err(9) = ',num2str(Error(9)),space,'|| s_min(9) = ',num2str(s_min(9))]);
        disp(['Err(10) = ',num2str(Error(10)),space,'|| s_min(10) = ',num2str(s_min(10))]);
        disp(['------------------------------------------------------------------------']);
        disp(' ');   
    end     
     
    %--- save ---%
     
    MIN_SV(:,:,p) = Min_sv;
    ERROR_BOX(:,:,p) = Error_box;
    N_VALUES(:,p) = N_values(p,:);
    M_VALUES(:,p) = M_values;
 
end       
 
%--- save result ---%
 
%save('Firstsam_d15_NlogN_r05','Error_box','N_values','M_values','Min_sv')
%save('Firstsam_d2to10_NlogN_r05','ERROR_BOX','N_VALUES','M_VALUES','MIN_SV')
t1 = toc;
 
%--- finish pararellize ---%
 
p = gcp;                                                              
delete(p)
delete(gcp('nocreate'))
