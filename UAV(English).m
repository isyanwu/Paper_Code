clc;
clear all;
delta_t=1;
% W=10^6;
W=10^3;
T=71;
tau=50;
N=T/delta_t;
M=tau/delta_t;
Fu_max=10^3;
Fg_max=10^5;
q_0=[125;125];
q_f=[125;125];
direction=([0;0]-q_0)/norm(([0;0]-q_0));
initialv=norm(q_f)/((T-1)/2);
v_0=initialv*direction;
v_f=-v_0;
V_max=30;
V_min=1;
A_max=12;
L=10^2;
H=10^2;
K=10^3;
P=0.1;
P=0.01;
beta_0=10^-6;
theta_2=10^-9;
gamma_0=beta_0*P/theta_2;
c1=9.26 *10^-4;
c2=2250;
% c2=1000;
g=9.8;
eta=10^-28;


%%生成初始解
q_j(:,1)=q_0;
v_j(:,1)=v_0;
for n=1:(N-1)/2
    v_j(:,n)=v_0;
    q_j(:,n+1)=q_j(:,n)+v_j(:,n)*delta_t;
end

for n=(N-1)/2+1:N-1
    v_j(:,n)=v_f;
    q_j(:,n+1)=q_j(:,n)+v_j(:,n)*delta_t;
end
v_j(:,N)=v_f;
% a=norm(2*v_j(:,1))


for n=1:N
omega_j(n)=V_min;
end
alphaa_j=0.5;
Fu_j=Fu_max;
Fg_j=Fg_max;
%

%验证初始解
a1=(1-alphaa_j)*K*L-Fu_j*T;
b=L*alphaa_j*K-Fg_j*(T-tau);
c=-1;
for n=1:N
    test=norm(v_j(:,n))-V_max;
    if test>0
        c=1;
    end
end

d=-1;
e=-1;
for n=1:N
    test=V_min-omega_j(n);
    if test>0
        d=1;
    end
    test=omega_j(n)^2-norm(v_j(:,n))^2;
    if test>0
        e=1;
    end
end

test=0;
for m=1:M
    test=test+W*log2(1+gamma_0/(H^2+norm(q_j(:,m))^2));
end
g1=alphaa_j*L-test;

h=-1;
i=-1;
for n=1:N-1    
    if v_j(:,n+1)~=v_j(:,n)+[0;0]*delta_t 
        h=1;%C5
    end
    if q_j(:,n+1)~=q_j(:,n)+v_j(:,n)*delta_t+[0;0]*delta_t^2/2;
        i=1;
    end%C6
end

%
a_j=zeros(2,length(v_j));
for i=1:length(v_j)-1
    a_j(:,i)=v_j(:,i+1)-v_j(:,i);
end
temp_1=0;
for n=1:N
%     c3=1/omega_j(n) - 1/(omega_j(n)^2)*(omega(n)-omega_j(n));
%     temp_1=temp_1+c1*pow_pos(norm(v(:,n)),3)+c2*c3 + 2*c2*pow_pos(norm(a(:,n)),2)/g^2/omega_j(n) - c2/(2*g^2*omega_j(n)^2)*(  (omega_j(n)+norm(a_j(:,n))^2)^2+2*(omega_j(n)+norm(a_j(:,n))^2)*(2*a_j(:,n)'*(a(:,n)-a_j(:,n))+omega(n)-omega_j(n)) - omega(n)^2 - pow_pos(norm(a(:,n)),4)     );
    temp_1=temp_1+c1*norm(v_j(:,n))^3+c2*1/omega_j(n)+c2/g^2*norm(a_j(:,n))^2/omega_j(n);
end
obj_j=temp_1+alphaa_j*eta*L*K*Fg_j^2+(1-alphaa_j)*eta*L*K*Fg_j^2+P*tau;
obj_j=obj_j/T;



%%cvx
f=[];
for k=1:50
cvx_begin quiet
cvx_precision best
variable alphaa
variable Fu
variable Fg
variable varphi
variable psii
variable q(2,N)
variable v(2,N)
variable a(2,N)
variable omega(N)
expression obj;
temp_1=0;
for n=1:N-1
%     c3=1/omega_j(n) - 1/(omega_j(n)^2)*(omega(n)-omega_j(n));
%     temp_1=temp_1+c1*pow_pos(norm(v(:,n)),3)+c2*c3 + 2*c2*pow_pos(norm(a(:,n)),2)/g^2/omega_j(n) - c2/(2*g^2*omega_j(n)^2)*(  (omega_j(n)+norm(a_j(:,n))^2)^2+2*(omega_j(n)+norm(a_j(:,n))^2)*(2*a_j(:,n)'*(a(:,n)-a_j(:,n))+omega(n)-omega_j(n)) - omega(n)^2 - pow_pos(norm(a(:,n)),4)     );
    temp_1=temp_1+c1*pow_pos(norm(v(:,n)),3)+c2*inv_pos(omega(n))+c2/g^2*quad_over_lin(a(:,n),omega(n));
end
obj=varphi+psii+temp_1;
% minimize varphi+psii+temp_1
minimize obj

subject to
0 < Fu <= Fu_max;       %C1
0< Fg <= Fg_max;       %C2
0<= alphaa <=1;         %C3
(1-alphaa)*K*L-Fu*T<=0;
L*alphaa*K-Fg*(T-tau)<=0;     %C4

for n=1:N-1    
    v(:,n+1)==v(:,n)+a(:,n)*delta_t;      %C5
    q(:,n+1)==q(:,n)+v(:,n)*delta_t+a(:,n)*delta_t^2/2;     %C6
end

q(:,1)==q_0;
q(:,N)==q_f;      %C7

v(:,1)==v_0;
v(:,N)==v_f;      %C8

for n=2:N-1
    norm(v(:,n))<=V_max;      %C9
end

for n=1:N-1
    norm(a(:,n))<=A_max;      %C10
end

for n=1:N
    V_min - omega(n)<=0;        %C11a
    omega(n)^2 - (  norm(v_j(:,n))^2 + 2*v_j(:,n)'*(v(:,n)-v_j(:,n))  )<=0;   %C11b
end

temp_2=0;
for m=1:M
    c3=H^2+norm(q_j(:,m))^2;
    temp_2=temp_2 + log2(1+gamma_0/c3) - gamma_0*(pow_pos(norm(q(:,m)),2) - norm(q_j(:,m))^2)/((c3+gamma_0)*c3*log(2));
end
alphaa*L/W- temp_2<=0;        %word里面写错了

1/2*eta*L*K*( pow_pos((alphaa+Fg^2),2) - (alphaa_j^2+2*alphaa_j*(alphaa-alphaa_j)) - (Fg_j^4+4*Fg_j^3*(Fg-Fg_j))  ) - varphi<=0;

1/2*eta*L*K*(   pow_pos(((1-alphaa)+Fu^2),2) - 1 + 2*alphaa -(alphaa_j^2+2*alphaa_j*(alphaa-alphaa_j)) - (Fu_j^4+4*Fu_j^3*(Fu-Fu_j))   ) - psii<=0;

cvx_end
q_j=q;
Fg_j=Fg;
Fu_j=Fu;
alphaa_j=alphaa;
v_j=v;

% temp_1=0;
% for n=1:N-1
% %     c3=1/omega_j(n) - 1/(omega_j(n)^2)*(omega(n)-omega_j(n));
% %     temp_1=temp_1+c1*pow_pos(norm(v(:,n)),3)+c2*c3 + 2*c2*pow_pos(norm(a(:,n)),2)/g^2/omega_j(n) - c2/(2*g^2*omega_j(n)^2)*(  (omega_j(n)+norm(a_j(:,n))^2)^2+2*(omega_j(n)+norm(a_j(:,n))^2)*(2*a_j(:,n)'*(a(:,n)-a_j(:,n))+omega(n)-omega_j(n)) - omega(n)^2 - pow_pos(norm(a(:,n)),4)     );
%     temp_1=temp_1+c1*pow_pos(norm(v(:,n)),3)+c2*inv_pos(omega(n))+c2/g^2*quad_over_lin(a(:,n),omega(n));
% end
% obj=varphi+psii+temp_1;
f = [f obj];
end
f=f+P*tau+1/2*m*(norm(v(:,N))^2-norm(v(:,1))^2);
figure(1);
plot(f,'-o')

figure(2);
x=q(1,1);
y=q(2,1);
for i=2:length(q)
    plot([x,q(1,i)],[y,q(2,i)]);
    hold on
    x=q(1,i);
    y=q(2,i);
end
% 
% %
% 
% test=0;
% for m=1:M
%     test=test+W*log2(1+gamma_0/(H^2+norm(q(:,m))^2));
% end
% error=alphaa*L-test;

save('..\model_1数据结果\p_alphaa','alphaa');
save('..\model_1数据结果\p_Fu','Fu');
save('..\model_1数据结果\p_Fg','Fg');
save('..\model_1数据结果\p_q','q');
save('..\model_1数据结果\p_v','v');
save('..\model_1数据结果\p_a','a');
save('..\model_1数据结果\p_f','f');


% save('..\model_1数据结果\h_alphaa','alphaa');
% save('..\model_1数据结果\h_Fu','Fu');
% save('..\model_1数据结果\h_Fg','Fg');
% save('..\model_1数据结果\h_q','q');
% save('..\model_1数据结果\h_v','v');
% save('..\model_1数据结果\h_a','a');
% save('..\model_1数据结果\h_f','f');
