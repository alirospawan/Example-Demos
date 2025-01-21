clear all;clc;

d=3;  % d+1=delay-time
T=1000; % T=epoch numbers
nu=d+2; % nu=input numbers
ny=3; % ny=output numbers
K=4; % K=node(rule) numbers
L=4; % L=enhancement node group numbers
m=4;
Np=7; % Np=predictive numbers;
M=3;
nL=M*K;%alpha=rand;

c=5*ones(M,K);
sigma=2*ones(M,K);


sp=[0*ones(d,1);0.1*ones(T/4,1);0.4*ones(T/4,1);0.5*ones(T/4,1);0.3*ones(T/4+M+d,1)];


y=0*ones(1,d+2);u=0*ones(1,d+2); 

R=0.002;
t = (0:T+10)';
E=R*randn(length(t),1);

KK=[0*ones(599,1); 0.1*ones(T+d-10,1)]; % KK:disturbance;

ayj(:,2)=0.6*ones(m,1);
p=0.08;
pp = 0.08;

error = 0;
sumError=0;
ym=0;

Wf=0.5*ones(K,1);
Wh=2.5*ones(m,1);     
Be=-0.257*ones(K,1);
Af=0.1*ones(K,1);
Ae=-2*ones(m,1);
We=0.05*ones(m,1);    % As small as possible


etaG=0.001;

q = 2;  % Elliptic shape
% For q=2, the membership function resembles a parabolic shape.
% For q>2, the membership function becomes flatter near the center and steeper near the edges.
% For q<2, the membership function becomes steeper near the center and flatter near the edges.

c_id=5*ones(11,11);
sigma_id=2*ones(11,11);

Wf_id=0.3*ones(11,1);
Wh_id=0.1*ones(11,1);     
Be_id=-1.5*ones(11,1);
Af_id=2*ones(11,1);
Ae_id=-1.5*ones(11,1);
We_id=0.01*ones(11,1);    % As small as possible


etaG_id = 0.001;
ym_id = 0;

tic;

for k=2:T
    error_old = error;
    error = sp(k)- y(k);
    sumError = error + sumError;

    ym_old = ym(k-1);

    x(1)= error; 
    x(2)= sumError;
    x(3)= error - error_old;

    mu=zeros(M,4);  
    tou_=ones(1,4); 
    for k0=1:K
        for t=1:M
            ratio = abs((x(t)-c(t,k0)) / sigma(t,k0)); 
            if ratio < 1
                mu(t,k0) = (1 - ratio^q)^(1 / q);  % Eq.(10)    % Elliptic TSK Fuzzy Membership Function
            else
                mu(t,k0) = 0.00001;
            end
        end
        for t=1:M
            tou_(k0) = tou_(k0) * mu(t,k0);
        end       
    end

    sum_tou_ = 0.0;
    for k0=1:K
        sum_tou_ = sum_tou_ + tou_(k0);
    end

    omega=zeros(1,4);
    for k0=1:K
        omega(k0)=tou_(k0)/sum_tou_;
    end

    alpha=zeros(M,4);
    for k0=1:K
        for t=1:M
            alpha(t,k0)= omega(k0);
        end
        % alpha(1,k0) = 0;
        % alpha(2,k0) = 0;
    end

    zk=zeros(3,4);
    for k0=1:K
        for t=1:M
            zk(t,k0)=alpha(t,k0)*x(t)+(Af(k0)*(p*ym_old)); 
        end
    end

    z_=zeros(1,4);  
    for k0=1:K
        z_(k0) = 0;
        for t=1:M
           z_(k0) = z_(k0) + zk(t,k0);  % z_i_sk
        end
        Z(k0)=omega(k0)*z_(k0);
    end

    omega_=[0 0 0 0];
    for k0=1:K
        omega_(k0) = omega_(k0) + omega(k0);
    end

    Fi=zeros(1,4);
    for k0=1:K
        Fi(k0)=z_(k0)*(omega_(k0))*Wf(k0);
    end

    F = 0;
    for k0=1:K
      F = F + Fi(k0);
    end
    

    xi=zeros(1,4);  H=zeros(1,4);
    for j=1:m
        for k0=1:K
            xi(j)=tanh(Z(k0)*We(j)+Be(j));
            H(j)= xi(j)*((Z(k0))*(We(j))+Be(j)+Ae(k0)*(p*ym_old));    %
        end
    end

    HWh = 0;
    for k0=1:K
      HWh =  HWh + (H(k0) * Wh(k0));
    end

    ym(k) = F + HWh ;

    ym_old =ym(k);
    


    % Adaptive Mechanism for System Identification START (id:Identification)
    ym_old_id = ym_id(k-1);

    for i = 1:ny+1
        if k - i >= 1
            x_idd(i) =  y(k-i+1);     % y(k), y(k-1), ..., y(k-n_y)
        else
            x_idd(i) = 0;
        end
    end
   
    for j = 1:nu+1
        if k - d - j >= 1
            x_idd(ny + j + 1) = ym(k - d - j + 1); % u(k-d), u(k-d-1), ..., u(k-d-n_u)
        else 
            x_idd(ny + j + 1) = 0;
        end
    end

    % y(k), y(k-1), y(k-2), y(k-3),  u(k-d), u(k-d-1),  u(k-d-2), u(k-d-3), u(k-d-4), u(k-d-5), v(k)
    x_id = [x_idd KK(k)];       % added disturbance to the model

    tou_id=ones(1,11); 
    mu_id=zeros(11,11);
    for k0=1:11
        for t=1:size(x_id,2)
            ratio_id = abs((x_id(t)-c_id(t,k0)) / sigma_id(t,k0)); 
            mu_id(t,k0) = (1 - ratio_id^q)^(1 / q);
        end
        for t=1:size(x_id,2)
            tou_id(k0) = tou_id(k0) * mu_id(t,k0);
        end       
    end

    sum_tou_id = 0.0;
    for k0=1:size(tou_id,2)
        sum_tou_id = sum_tou_id + tou_id(k0);
    end

    omega_id=zeros(1, size(tou_id,2));
    for k0=1:size(tou_id,2)
        omega_id(k0)=tou_id(k0)/sum_tou_id;
    end

    alpha_id=zeros(size(x_id,2), size(tou_id,2));
    for k0=1:K
        for t=1:size(x_id,2)
            alpha_id(t,k0)= omega_id(k0);
        end
    end

    zk_id=zeros(size(x_id,2),size(tou_id,2));
    for k0=1:K
        for t=1:size(x_id,2)
            zk_id(t,k0)=alpha_id(t,k0)*x_id(t)+(Af_id(k0)*(pp*ym_old_id)); 
        end
    end

    z_id=zeros(1, size(tou_id,2));  
    for k0=1:size(tou_id,2)
        z_id(k0) = 0;
        for t=1:size(x_id,2)
           z_id(k0) = z_id(k0) + zk_id(t,k0);  % z_i_sk
        end
        Z_id(k0)=omega_id(k0)*z_id(k0);
    end
    omega__id=[0 0 0 0 0 0 0 0 0 0 0];
    for k0=1:size(tou_id,2)
        omega__id(k0) = omega__id(k0) + omega_id(k0);
    end

    Fi_id=zeros(1, size(tou_id,2));
    for k0=1:size(tou_id,2)
        Fi_id(k0)=z_id(k0)*(omega__id(k0))*Wf_id(k0);
    end

    F_id = 0;
    for k0=1:size(tou_id,2)
      F_id = F_id + Fi_id(k0);
    end

    xi_id=zeros(1, size(tou_id,2));  H_id=zeros(1, size(tou_id,2));
    for j=1:m
        for k0=1:size(tou_id,2)
            xi_id(j)=tanh(Z_id(k0)*We_id(j)+Be_id(j));
            H_id(j)= xi_id(j)*((Z_id(k0))*(We_id(j))+Be_id(j)+Ae_id(k0)*(pp*ym_old_id));    %
        end
    end

    HWh_id = 0;
    for k0=1:size(tou_id,2)
      HWh_id =  HWh_id + (H_id(k0) * Wh_id(k0));
    end

    ym_id(k) = F_id + HWh_id ;

    %|% Error Objective function = 0.5 * ym_id(k)-y(k) ^2;

    dE_dWf_id = (ym_id(k)-y(k)).*Z_id';
    Wf_id   = Wf_id - etaG_id *dE_dWf_id;

    dE_dWh_id = (ym_id(k)-y(k)).*H_id';
    Wh_id   = Wh_id - etaG_id *dE_dWh_id;

    dE_dBe_id = (ym_id(k)-y(k)).*Wh_id.*(1-H_id.*H_id)';
    Be_id   = Be_id - etaG_id * dE_dBe_id;

    dE_dAf_id = (ym_id(k)-y(k)).*Wf_id.*ym(k-1).*(1-Z_id.*Z_id)';
    Af_id   = Af_id - etaG_id * dE_dAf_id;

    dE_dAe_id = (ym_id(k)-y(k)).*Wh_id.*ym(k-1).*(1-H_id.*H_id)';
    Ae_id   = Ae_id - etaG_id * dE_dAe_id;

    dE_dWe_id = (ym_id(k)-y(k)).*Wh_id.*Z_id'.*(1-H_id.*H_id)';
    We_id   = We_id - etaG_id * dE_dWe_id;

    cont_weights_id(:,k) = [Wf_id(1); Wh_id(1); Be_id(1); Af_id(1); Ae_id(1); We_id(1)];

    norm_dydWf_id = norm(Z_id)^2;        % norm_X_squared_manual = sum(X.^2);
    norm_dydWh_id = norm(H_id)^2;
    norm_dydBE_id = norm(Wh_id.*(1-H_id.*H_id))^2;
    norm_dydAf_id = norm(Wf_id.*ym(k-1).*(1-Z_id.*Z_id)')^2;
    norm_dydAe_id = norm(Wh_id.*ym(k-1).*(1-H_id.*H_id)')^2;
    norm_dydWe_id = norm(Wh_id.*Z_id'.*(1-H_id.*H_id)')^2;

    max_norm_id = max([norm_dydWf_id norm_dydWh_id norm_dydBE_id norm_dydAf_id norm_dydAe_id norm_dydWe_id]);

    beta_id = 1;
    etaG_id = beta_id / (max_norm_id);
    eta_id(k+1) = beta_id / (max_norm_id);

    % Adaptive Mechanism for System Identification END

    e_cost = sp(k) - ym_id(k);

    dE_dWf = e_cost.*Z';
    Wf   = Wf - etaG *dE_dWf;

    dE_dWh = e_cost.*H';
    Wh   = Wh - etaG *dE_dWh;

    dE_dBe = e_cost.*Wh.*(1-H.*H)';
    Be   = Be - etaG * dE_dBe;

    dE_dAf = e_cost.*Wf.*ym(k-1).*(1-Z.*Z)';
    Af   = Af - etaG * dE_dAf;

    dE_dAe = e_cost.*Wh.*ym(k-1).*(1-H.*H)';
    Ae   = Ae - etaG * dE_dAe;

    dE_dWe = e_cost.*Wh.*Z'.*(1-H.*H)';
    We   = We - etaG * dE_dWe;

    cont_weights(:,k) = [Wf(1); Wh(1); Be(1); Af(1); Ae(1); We(1)];

    norm_dydWf = norm(Z)^2;        % norm_X_squared_manual = sum(X.^2);
    norm_dydWh = norm(H)^2;
    norm_dydBE = norm(Wh.*(1-H.*H))^2;
    norm_dydAf = norm(Wf.*ym(k-1).*(1-Z.*Z)')^2;
    norm_dydAe = norm(Wh.*ym(k-1).*(1-H.*H)')^2;
    norm_dydWe = norm(Wh.*Z'.*(1-H.*H)')^2;

    max_norm = max([norm_dydWf norm_dydWh norm_dydBE norm_dydAf norm_dydAe norm_dydWe]);



    beta = 0.5;
    etaG = beta / (max_norm);
    eta(k+1) = beta / (max_norm);

    delta_u(k)=ym(k)-ym(k-1);
    u(k) =  u(k-1) + delta_u(k);

    ud=[0*ones(1,d-1) u]; ud1=[0*ones(1,d) u];
    y(k+1)=y(k)^3-0.2*abs(y(k))*ud(k)+0.08*ud(k)^2+E(k)+KK(k);

end
toc;

% find the performance index
e_y=abs(sp(1:T)'-y(1:T));
maxE_y=max(e_y);
rmse_y=sqrt((e_y*e_y')/T);
ISE_y=sum(e_y*e_y');
IAE_y=sum(e_y);

for ttt=1:T
    xx(ttt)=ttt*e_y(ttt);
end
ITAE_y=sum(xx);


fprintf('maxError: %.4f\n', maxE_y);
fprintf('RMSE: %.4f\n', rmse_y);
fprintf('ISE: %.4f\n', ISE_y);
fprintf('IAE: %.4f\n', IAE_y);
fprintf('ITAE: %.4f\n', ITAE_y);


tvec=1:T;

figure(1);
plot(tvec,y(1:T),'r',tvec,sp(1:T),'b:')
xlabel('sampling number');
ylabel('sp & y');
legend('PLANT OUTPUT','SETPOINT','Location','southeast')


figure(2);
subplot(211)
plot(tvec,ym_id(1:T),'r');xlabel('sampling number');ylabel('y identifier');
subplot(212)
plot(tvec,ym(1:T),'r');xlabel('sampling number');ylabel('y controller');


figure(3);
plot(tvec,ym_id(1:T),'b', tvec,ym(1:T),'r'); xlabel('sampling number');ylabel('y');
legend('y identifier','y controller')

