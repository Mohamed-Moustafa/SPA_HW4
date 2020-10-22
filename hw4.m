%% **************************** part 1 ************************* %%
clear
close all
% ************************** 1 - givens *************************** %

% A matrix
A=[0.7 0.5 0;
   -0.5 0.7 0;
   0 0 0.9];

% B vector
B=[1;1;1];

% C vector
C=[0 -1 1];

%D=0.5
D=0.5;

% N
N=0:19;

% y ref
y_ref=[1 0 0 4 4 1 0 2 0 1 3 4 4 2 1 2 4 3 2 2]';

% xo
xo=[0.1;0.2;0.3];



%% ******************** 2 - Calculate Q *************************** %

% let's compute the Qk and put them in a vector
% Q_vector= [Qo Q1 Q2 .... Qn-1]
Q_vector=zeros(1,20);

% Qo = 0.5
Q_vector(1,1)=D;


for i=2:20
    Q_vector(1,i)= C* A^(i-2) * B;
end


% let's estaplish our Q matrix

Q=zeros(20,20);
for ii=1:20
    L=length(Q_vector)-(ii-1);
    vec= Q_vector(1,ii)*ones(1,L);
    
    Q_k= diag(vec,-(ii-1));
    Q= Q+Q_k;
end

%% **************** 3 - calculate phi & u_opt ******************* %
phi=zeros(20,3);
for i=1:20
   phi(i,:)= C * A^(i-1) ;

end

u_opt = (Q' * Q) \ Q' * (y_ref - phi*xo);

% here we limit our control output from -100 to 100
u_limit=u_opt;
u_limit(u_limit>100)=100;
u_limit(u_limit<-100)=-100;



%% **************** 4 - calculate xk & yk ******************* %


xk_2=zeros(3,20);
yk_2=zeros(1,20);

xk_3=zeros(3,20);
yk_3=zeros(1,20);

xk_2(:,1)= xo;
for i=1:20
   %phi(i,:)= C * A^(i-1) ;
    %xk(:,i)=  A^(i-1)*xo ;%+ B*u_opt(i) ;
    if i~=1
        xk_2(:,i)= A*xk_2(:,i-1)+B*u_opt(i-1);
        
    end
    yk_2(i)= C * xk_2(:,i)+D*u_opt(i);
end



%% **************** 5 - calculate yk & yk_limit ******************* %

yk= phi * xo + Q * u_opt;
yk_limit= phi * xo + Q * u_limit;


%% **************** 6 - plot ****************************** %
% 
% plot yk and y_ref on same graph
figure;hold on;
plot(N,yk,'b');
plot(N,y_ref,'g--');
title('Comparing yK with yRef (u* case) ');
legend('yK','yRef')
xlabel('Horizon (N)','FontSize',12);
ylabel('Trajectory ref','FontSize',12);

% plot yk_limit and y_ref 
figure;hold on;
plot(N,yk_limit);
plot(N,y_ref,'r')
title('Comparing yK with yRef (uLimit case) ');
legend('yK','yRef')
xlabel('Horizon (N)','FontSize',12);
ylabel('Trajectory ref','FontSize',12);
grid on

%plot the control outputs u_optimal and u_limit
figure;hold on;
plot(N,u_opt,'r');
plot(N,u_limit ,'b')
title('our control signal');
legend('u*','uLimit')
xlabel('Horizon (N)','FontSize',12);
ylabel('Control signal','FontSize',12);
grid on

%% **************************** part 2 ************************* %%


% 1- Design our matrcies Phi , q , H , R , P

% phi matrix (Transition State Matrix)
phi=[0.7 0.5 0;
   -0.5 0.7 0;
   0 0 0.9];

% Q matrix (White Noise Covariance matrix)
q=   [0.5 0 0; 
      0 0.5 0;
      0 0 0.5];
% H matrix
H=   [1 0 0; 
      0 1 0;
      0 0 1];
% R matrix (Sensor Covariance matrix)
R=  [1.7 0 0; 
     0 1 0;
     0 0 1.8];
  
% P matrix
P=   [1 0 0; 
      0 1 0;
      0 0 1];
  
  
  
%% 2- now let's get our zk (measurments data) which are zk = xk + [N (0, 1.7), N (0, 1.0), N (0, 1.8)]T

zk=zeros(3,20);
for i=1:20
zk(:,i)= xk_2(:,i) + [1.7*randn(); randn() ; 1.8*randn()];
end


%% 3- now let's start our Kalman filter loop

%inital value for X vector
X(:,1)=xk_2(:,1);

%the obtained y (output) of kalman intalization
y_kalman=zeros(1,20);

%the kalman_gain which are of size 3x3x20 ---> 3x3 is the gain matrix and
%20 for horizion
Kalman_gain=zeros(3,3,20);

%now intializing the first element of by 3x3 zeros matrix
Kalman_gain(:,:,1)=zeros(3,3);



for i=2:20
   
   Zk=zk(:,i);
   
   %prediction step
   Xp= phi* X(:,i-1)+B*u_opt(i-1);
   Pp= phi * P * phi' + q;
   
   %update step
   S= H * Pp * H' + R;
   K= Pp * H' * inv(S);
   Y(:,i)= Zk - H * Xp;
   
   X(:,i)= Xp + K * Y(:,i);
   P= (eye(3) - K*H) * Pp * (eye(3) - K*H)' + K* R *K';
   Kalman_gain(:,:,i)= K;
   
end

% Calculating the y_obtained from kalman
for j=1:20
    y_kalman(1,j)=C*X(:,j)+D*u_opt(j);
 
end

%% 5- now let's plot our graphs




figure; plot(N,y_kalman,'b');

hold on
plot(N,y_ref,'g--');

title('yKalman & yRef');
legend('yKalman','yRef')
xlabel('Horizon (N)','FontSize',12);
ylabel('Trajectory','FontSize',12);
grid on

figure; plot(N, squeeze(Kalman_gain(1,1,:)))
hold on
plot(N,squeeze(Kalman_gain(2,2,:)));
plot(N,squeeze(Kalman_gain(3,3,:)));
title('Kalman gain');
legend('k(1,1)','k(2,2)','k(3,3)')
xlabel('Horizon (N)','FontSize',12);
ylabel('gain','FontSize',12);
grid on
