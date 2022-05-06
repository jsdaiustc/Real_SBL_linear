% Matlab code for "Real-Valued Sparse Bayesian Learning for DOA Estimation With Arbitrary Linear Arrays"

function [Pm,search_area]=real_VBI(X,K,resolution,position,etc)
search_area=[-90:resolution:90];

%% real-valued transformation
Rx=X*X'/K;
Y=Rx(:);
[M2,~]=size(Y);
M=sqrt(M2);
K_hat=length(search_area);
In=eye(M);
In=In(:);
pos_all=round(log(kron(exp(-position),exp(position))));
W=kron(Rx.',Rx)/K;
W_sq=sqrtm(inv(W));
Y=W_sq*Y;
Y=real(Y);
a_search=[-90:2:90]*pi/180;
A=exp(-1i*pi*pos_all'*sin(a_search))/sqrt(1);
A_w=W_sq*A;
[Ve,De]=eig(real(A_w'*A_w));
de=diag(De);
ind_d1=find(de > de(end)/1e8);
Ve=Ve(:,ind_d1);
U1= A_w*Ve*  diag(de(ind_d1).^(-1/2));
W_sq=U1'*W_sq;
Y=real(U1'*Y);
M_eff=length(ind_d1);

%% Initialization for VBI
a=1e-5;
b=1e-5;
maxiter=300;
tol=1e-5;
beta0=.01;
delta_inv=ones(K_hat+1,1);
converged = false;
iter = 0;
a_search=search_area*pi/180;
A=exp(-1i*pi*pos_all'*sin(a_search));
B=-1i*pi*pos_all'*cos(a_search).*A;
A_w=real(W_sq*A);
B_w=real(W_sq*B);
I_w= real(W_sq* In);
Phi=[A_w, I_w  ];


%% Bayesian Inference
while (~converged) || iter<=100
    
    iter = iter + 1;
    delta_last = delta_inv;

    %% calculate mu and Sigma
    V_temp=  1/beta0*eye(M_eff) +  Phi*diag(delta_inv) * Phi';
    Sigma = diag(delta_inv) -diag(delta_inv) * Phi' * (V_temp\Phi) * diag(delta_inv);
    mu = beta0*(Sigma*(Phi'*Y));
     
    %% update delta_inv
    bb= abs(mu).^2 +    diag(Sigma);
    c_k = a + 1/2;
    d_k = b + bb/2;
    delta_inv=d_k./c_k;
    
    %% stopping criteria
    erro=norm(delta_inv - delta_last)/norm(delta_last);
    if erro < tol || iter >= maxiter
        converged = true;
    end
    
    %% grid refinement
    sum_mu=sum( mu.*conj(mu), 2);
    Pm=sum_mu(1:end-1);
    [~, idx] = sort(Pm(1:end-1), 'descend');
    idx = idx(1:etc);
    BHB = B_w' * B_w;
    P = real(conj(BHB(idx,idx)) .* (mu(idx,:) * mu(idx,:)' + Sigma(idx,idx)));
    v =  real(diag(conj(mu(idx))) * B_w(:,idx)' * (Y - A_w * mu(1:end-1)-mu(end)*I_w))...
        -  real(diag(B_w(:,idx)' * A_w * Sigma(1:end-1,idx))  +    diag(Sigma(idx,K_hat+1))*B_w(:,idx).'*conj(I_w));
    temp1=v./diag(P);
    temp2=temp1'*180/pi;
    if iter<100
        ind_small=find(abs(temp2)<resolution/100);
        temp2(ind_small)=sign(temp2(ind_small))*resolution/100;
    end
    ind_large=find(abs(temp2)>resolution);
    temp2(ind_large)=sign(temp2(ind_large))*resolution/100;
    angle_cand=search_area(idx) + temp2;
    search_area(idx)=angle_cand;
    A_ect=exp(-1i*pi*pos_all'*sin(search_area(idx)*pi/180));
    B_ect=-1i*pi*pos_all'*cos(search_area(idx)*pi/180).*A_ect;
    A_w(:,idx) =real(W_sq*A_ect);
    B_w(:,idx) =real(W_sq*B_ect);
    Phi(:,idx)= A_w(:,idx);
    
    if iter==50
        beta0=1;
    end
end

%% Output
Pm=delta_inv(1:end-1);

%[search_area,sort_s]=sort(search_area);
% Pm=Pm(sort_s);
% insert=(search_area(1:end-1)+search_area(2:end))/2;
% search_area_2=zeros(length(search_area)*2-1,1);
% search_area_2(1:2:end)=search_area;
% search_area_2(2:2:end)=insert;
% Pm_2=zeros(length(Pm)*2-1,1);
% Pm_2(1:2:end)=Pm;
% search_area=search_area_2;
% Pm=Pm_2;
