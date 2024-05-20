%-------------------------------------------------------------------------%
% Optimal Sampling Version 3
% Description: We solve the weighted-LSP approximation of a function on an
% irregular domain, from sequential samples point.
% Programer: Juan Manuel Cardenas
% Date: July 9 - 2019 / Last modification: July 16 - 2019
%-------------------------------------------------------------------------%
 
%% Set up
tic
clear all; close all; clc;
space=' ';
%parpool;                        % Start parallelize 
index_type = 1;                 % Index set HC
d_values =  [15]; % [2 3 5 10]; % Dimension
m = 0;                          % number of k_points: M = m*N or M=N*log(N)
domain = -1;                    % 0 cube\corner, 1 cube\sphere, -1 annular  
f = 1;                          % type of function: f=1 , f= 0 Genz
Trials = 10;                    % number of trials
rad = 1/2;                        % annular's radio
N_max= 1000;                    % biggest size of |HC|
max_div = 6;                   % nmax with d=15 is 6, with d = 10 and 10>d is 20.
 
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
         
            x = 2*rand(1,d)-rand(1,d);
         
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
 
    %--- Measure mu---%
     
    mu = abs(Q.^2);                               
     
    %--- Function right-hand side ---%
     
    if f == 1
       f_grid = func(Z);
    else
       f_grid = F_Genz(Z);
    end
     
%%  Compute W-LSP approximation      
 
    %--- M points ---%
     
    if m == 0                          
       k_pts = round(log(N_values(p,:)));                    % M = kN 
       M_values = k_pts.*N_values(p,:);
    else
       k_pts = m*N;                    
       M_values = k_pts.*N_values(p,:);
    end
     
    s_min = [];
    Error = [];
    I_new = zeros(max(k_pts),max(N_values(p,:)));
    I_ad = [];
     
    %--- Trials ---%
     
    for t = 1 : Trials
         
        AD = 0;
         
        for l = 1:length(N_values(p,:))
             
            N = N_values(p,l);
            M = M_values(l);
             
            %--- Choosing M points i.i.d under d\mu ---%
             
            for j = 1 + AD : N
                I_new(1:k_pts(l),j) = datasample(1:K,k_pts(l),'Replace',true,'Weights',mu(:,j));
            end
             
            if  ( 1  < l) && ( l < length(k_pts) ) 
                k_ad = k_pts(l) - k_pts(l-1);
                I_ad = zeros(k_ad,N_values(p,l-1));
                for j = 1: N_values(p,l-1)
                    I_ad(1:k_ad,j) = datasample(1:K,k_ad,'Replace',true,'Weights',mu(:,j));
                end  
            end
             
            I_NEW = nonzeros( I_new(:,1+AD:N) );
            I_AD = nonzeros( I_ad );
            AD = N;
             
            if N == N_values(p,1)
                I = I_NEW;
            else
                I = [I ; I_AD ; I_NEW];
            end 
             
            %--- Data matrix ---%
 
            W = (1/N)*(sum((mu(I,1:N))')).^(-1);   % Compute the weights  
            W_D = diag(sqrt(W));                   % Save weights associate to M pts 
            A = (sqrt(K)/sqrt(M))*W_D*Q(I,1:N); 
             
            %--- Right-hand side ---%
         
            F = f_grid(I);
            b = (1/sqrt(M))*W_D*F;
 
            %--- Compute least-squares approximation ---%
         
            c = A\b;
            s_min(l) = min(svd(A));                    % save min singular value
         
            %--- Compute error ---%
         
            f_tilde = sqrt(K)*Q(:,1:N)*c;
            Error(l) = norm(f_grid - f_tilde);
         
            disp(' ');
            disp(['--------------------------------- DONE ---------------------------------']);
            disp(['dimension = ',num2str(d),space,'|| N = ',num2str(N),space,'|| Iteration = ',num2str(l),space,'|| Measure = Adaptative',space,'||Trial = ',num2str(t)]);
            disp(['Error = ',num2str(Error(l)),space,'|| s_min = ',num2str(s_min(l))]);
            disp(['------------------------------------------------------------------------']);
            disp(' ');
        end
         
        %--- save ---%
         
        Min_sv(:,t) = s_min;
        Error_box(:,t) = Error;
         
    end 
     
    %--- save ---%
     
    MIN_SV(:,:,p) = Min_sv;
    ERROR_BOX(:,:,p) = Error_box;
    N_VALUES(:,p) = N_values(p,:);
    M_VALUES(:,p) = M_values;
     
end
 
%--- save result ---%
 
%save('Newsam_d15_NlogN_r05','ERROR_BOX','N_VALUES','M_VALUES','MIN_SV')
%save('Newsam_d2to10_NlogN_r05','ERROR_BOX','N_VALUES','M_VALUES','MIN_SV')
t1 = toc;
 
%--- finish pararellize ---%
% p = gcp;                                                              
% delete(p)
% delete(gcp('nocreate'))
