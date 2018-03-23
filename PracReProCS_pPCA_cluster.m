% Han Guo, Chenlu Qiu, Namrata Vaswani 
% Copyright: Iowa State University
% Questions: hanguo@iastate.edu
% Reference: An Online Algorithm for Separating Sparse and Low-dimensional Signal Sequences from their Sum,
%            accepted by IEEE Transaction on Signal Processing


clear all; clc;
addpath Yall1;
addpath Data;
addpath Utils
addpath CS
addpath cluster_approx
load DarkVideoData_small.mat
DataTrain=Data(:,1:526);
I=Data(:,527:626);
Kmin     = 3;
Kmax     = 10;
alpha    = 20;
b        = 0.95;
p = size(I,1);
q = size(I,2);
imSize=imgsize;

%%%% training
mu0         = mean(DataTrain,2);
numTrain    = size(DataTrain,2);
MTrain      = DataTrain-repmat(mu0,1,numTrain);
M           = I-repmat(mu0,1, size(I,2)); %subtract the mean
[U, Sig, ~] = svd(1/sqrt(numTrain)*MTrain,0);

evals       = diag(Sig).^2;
energy      = sum(evals);
cum_evals   = cumsum(evals);
ind0        = find(cum_evals < b*energy);
rhat        = min(length(ind0),round(numTrain/10));
lambda_min  = evals(rhat);
U0          = U(:, 1:rhat); 


%%%% practical-ReProCS

Shat_mod    = zeros(p,q); 
Lhat_mod    = zeros(p,q); 
Nhat_mod    = cell(q,1); 
Fg          = zeros(p,q);
xcheck      = zeros(p,q);
Tpred       = [];
Ut          = U0; 

Pstar    = Ut;
k        = 0;
K        = [];
addition = 0;
cnew     = [];
t_new    = []; time_upd = []; thresh_diff = []; thresh_base = [];

for t = 1: 100
  t
clear opts; 

opts.selectAtom = 2;
opts.gamma = 0.1;
opts.tol=1e-4;
opts.maxIter=50;
opts.verbose=0;
opts.tol_early=0.005;
opts.debias = 0;
opts.hhs = 0;
opts.imheight = imSize(1);
opts.imwidth = imSize(2);
opts.display_perf = 1;
opts.SF = 0.25;
opts.thr = 0.05;
opts.iter = 20;
opts.display_progress = 0;
recflag = 'cluster';

Atf.times  = @(x) Projection(Ut,x); Atf.trans = @(y) Projection(Ut,y);
yt         = Atf.times(M(:,t));
k=1120;
c=1;
[xp,tcgs,tmodel,tproj] = cosamp_cluster(recflag, yt, k, c, Atf.times, Atf.trans, opts);




    
    omega(t) = sqrt(M(:,t)'*M(:,t)/p);
    That=find(abs(xp)>=omega(t));
    Shat_mod(That,t)    = subLS(Ut,That,yt);     
    Lhat_mod(:,t)       = M(:,t) - Shat_mod(:,t);
    Fg(That,t)          = 1;
    Nhat_mod{t}         = That;
    Tpred               = That;
    Bg(:,t)             = Lhat_mod(:,t) + mu0;

%     %% Projection PCA   
%     if addition==0    %&& norm( (Lhat(:,t-alpha+1:t) - Phat*(Phat'*Lhat(:,t-alpha+1:t)))./sqrt(alpha) )>thresh
%         addition = 1;
%         t_new    = t;
%         Pnewhat  = [];
%         k        = 0;
%     end
%         
%     if addition==1&&mod(t-t_new+1,alpha)==0
%         time_upd = [time_upd,t];           
%         D        = Lhat_mod(:,t-alpha+1:t)-Pstar*(Pstar'*Lhat_mod(:,t-alpha+1:t)); 
%        
%         [Pnew_hat, Lambda_new,~] = svd(D./sqrt(alpha),0);
%         Lambda_new               = diag(Lambda_new).^2;
%         Lambda_new               = Lambda_new(Lambda_new>=lambda_min);
%         th                       = round(rhat/3);
%         if size(Lambda_new,1)> th
%             Lambda_new=Lambda_new(1:th);
%         end
%            if numel(Lambda_new)==0
%                addition  = 0; 
%                cnew      = [cnew 0];
%            else              
%                cnew_hat    = numel(Lambda_new);
%                Pnewhat_old = Pnewhat;
%                Pnewhat     = Pnew_hat(:,1:cnew_hat); cnew = [cnew cnew_hat];%size(Pnewhat,2)];
%                Ut          = [Pstar Pnewhat];   
%                
%                k=k+1;
%                
%                if k==1 
%                    temp        =(Pnewhat*(Pnewhat'*Lhat_mod(:,t-alpha+1:t)));
%                    thresh_base = [thresh_base norm(temp./sqrt(alpha))];
%                    thresh_diff = [thresh_diff norm(temp./sqrt(alpha))];                 
%                else
%                    temp        =(Pnewhat*(Pnewhat'*Lhat_mod(:,t-alpha+1:t)));
%                    thresh_base = [thresh_base norm(temp./sqrt(alpha))];
%                    
%                    temp        = (Pnewhat*(Pnewhat'*Lhat_mod(:,t-alpha+1:t)) - Pnewhat_old*(Pnewhat_old'*Lhat_mod(:,t-alpha+1:t)));
%                    thresh_diff = [thresh_diff norm(temp./sqrt(alpha))];  
%                end
%                
%                flag = 0;
%                if k >= Kmin
%                    numK = 3;
%                    flag = thresh_diff(end)/thresh_base(end-1)<0.01;
%                    for ik = 1:numK-1
%                        flag = flag && thresh_diff(end-ik)/thresh_base(end-ik-1)<0.01;
%                    end
%                end
%                
%                if  k==Kmax|| (k>=Kmin && flag==1)                  
%                    addition =0;
%                    K        = [K k];
%                    Pstar    = Ut;            
%                end
%            end
%     end
%  
end
   
   

