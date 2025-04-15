format long
f=@(t,y) exp(-y)*sin(t+(2*pi*y));
%the ODE we will be evaluating
T=0;
Y=1;
%initial values
q(1:11,1)=1;
w(1:11,1)=1;
e(1:11,1)=1;
r(1:11,1)=1;
%these are matrices that will store values for each method
%this initial value is for Y(0)=1
for j=2:12
%repeating the methods for each n
n=2^(j);
h=2/n;
for i=1:n
    %evaluation of RK2
    K1=h*f(T,Y);
    K2=h*f(T+h,Y+K1);
    Y=Y+(K1/2)+(K2/2);
    %our Y for each step
    T=T+h;
    %our T for this step
    q(j-1,i+1)=Y;
    %the matrix that stores the Y-coordinates
end
RK2(j-1)=Y;
%the solution at T=2 for each n value
T=0;
Y=1;
for i=1:n
    K1=h*f(T,Y);
    K2=h*f(T+h,Y+K1);
    Y=Y+(0.45*K1)+(0.55*K2);
    T=T+h;
    w(j-1,i+1)=Y;
    %these are the same steps as RK2 above, except with B1=0.45 and B2=0.55
end
MRK2(j-1)=Y;
%the solution at T=2 for each n value
T=0;
Y=1;
for i=1:n
    Y=Y+(h*f(T,Y));
    T=T+h;
    e(j-1,i+1)=Y;
    %evaluation of Euler's method
end
EM(j-1)=Y;
%the solution at T=2 for each n value
T=0;
Y=1;
for i=1:n
    K1=h*f(T,Y);
    K2=h*f(T+(h/2),Y+(K1/2));
    K3=h*f(T+(h/2),Y+(K2/2));
    K4=h*f(T+h,Y+K3);
    Y=Y+((K1+(2*(K2+K3))+K4)/6);
    T=T+h;
    r(j-1,i+1)=Y;
    %evaluation of RK4
end
RK4(j-1)=Y;
%the solution at T=2 for each n value
T=0;
Y=1;
end
h=2/8192;
%h for evaluating the true solution for RK4 where n=8192
for i=1:8192
K1=h*f(T,Y);
K2=h*f(T+(h/2),Y+(K1/2));
K3=h*f(T+(h/2),Y+(K2/2));
K4=h*f(T+h,Y+K3);
Y=Y+((K1+(2*(K2+K3))+K4)/6);
T=T+h;
%evaluation of RK4 at n=8192
end
A=Y;
%our true solution
for i=1:11
ERK2(i)=log10(abs(RK2(i)-A));
EMRK2(i)=log10(abs(MRK2(i)-A));
EEM(i)=log10(abs(EM(i)-A));
ERK4(i)=log10(abs(RK4(i)-A));
%calculating the log of error at T=2 for each n and each method
end
for i=1:11
LOGN(i)=log10(2^(i+1));
%the log of n number of steps
end
hold on
figure(1)
plot(LOGN,ERK2,'r')
plot(LOGN,EMRK2,'b')
plot(LOGN,EEM,'g')
plot(LOGN,ERK4,'y')
%plotting the logs of the errors
hold off
for i=[2:12]
   figure(i)
   hold on
   plot([0:2/(2^i):2],q(i-1,[1:(1+(2^i))]),'r')
   plot([0:2/(2^i):2],w(i-1,[1:(1+(2^i))]),'b')
   plot([0:2/(2^i):2],e(i-1,[1:(1+(2^i))]),'g')
   plot([0:2/(2^i):2],r(i-1,[1:(1+(2^i))]),'y')
   %plotting the solutions to each method. Each n number of steps is a new graph
   hold off
end
%below we are just labeling our graphs
figure(1)
legend('RK2 error','My own RK2 error','EM error','RK4 error')
lgd=legend;
lgd.Title.String = 'Log of Errors';
figure(2)
legend('RK2','My RK2','EM','RK4')
lgd=legend;
lgd.Title.String = 'n=4 steps';
figure(3)
legend('RK2','My RK2','EM','RK4')
lgd=legend;
lgd.Title.String = 'n=8 steps';
figure(4)
legend('RK2','My RK2','EM','RK4')
lgd=legend;
lgd.Title.String = 'n=16 steps';
figure(5)
legend('RK2','My RK2','EM','RK4')
lgd=legend;
lgd.Title.String = 'n=32 steps';
figure(6)
legend('RK2','My RK2','EM','RK4')
lgd=legend;
lgd.Title.String = 'n=64 steps';
figure(7)
legend('RK2','My RK2','EM','RK4')
lgd=legend;
lgd.Title.String = 'n=128 steps';
figure(8)
legend('RK2','My RK2','EM','RK4')
lgd=legend;
lgd.Title.String = 'n=256 steps';
figure(9)
legend('RK2','My RK2','EM','RK4')
lgd=legend;
lgd.Title.String = 'n=512 steps';
figure(10)
legend('RK2','My RK2','EM','RK4')
lgd=legend;
lgd.Title.String = 'n=1024 steps';
figure(11)
legend('RK2','My RK2','EM','RK4')
lgd=legend;
lgd.Title.String = 'n=2048 steps';
figure(12)
legend('RK2','My RK2','EM','RK4')
lgd=legend;
lgd.Title.String = 'n=4096 steps';


format long
f=@(t,y) exp(-y)*sin(t+(2*pi*y));
%the ODE we will be evaluating
T=0;
Y=1;
%initial values
for j=2:12
%repeating the methods for each n
n=2^(j);
h=2/n;
for i=1:n
    %evaluation of RK2
    K1=h*f(T,Y);
    K2=h*f(T+h,Y+K1);
    Y=Y+(K2/2)+(K1/2);
    %our Y for each step
    T=T+h;
    %our T for this step
end
RK2(j-1)=Y;
%the solution at T=2 for each n value
T=0;
Y=1;
for i=1:n
    K1=h*f(T,Y);
    K2=h*f(T+h,Y+K1);
    Y=Y+(0.45*K1)+(0.55*K2);
    T=T+h;
    %these are the same steps as RK2 above, except with B1=0.45 and B2=0.55
end
MRK2(j-1)=Y;
%the solution at T=2 for each n value
T=0;
Y=1;
for i=1:n
    Y=Y+(h*f(T,Y));
    T=T+h;
    %evaluation of Euler's method
end
EM(j-1)=Y;
%the solution at T=2 for each n value
T=0;
Y=1;
for i=1:n
    K1=h*f(T,Y);
    K2=h*f(T+(h/2),Y+(K1/2));
    K3=h*f(T+(h/2),Y+(K2/2));
    K4=h*f(T+h,Y+K3);
    Y=Y+((K1+(2*(K2+K3))+K4)/6);
    T=T+h;
    %evaluation of RK4
end
RK4(j-1)=Y;
%the solution at T=2 for each n value
T=0;
Y=1;
end
h=2/8192;
%h for evaluating the true solution for RK4 where n=8192
for i=1:8192
K1=h*f(T,Y);
K2=h*f(T+(h/2),Y+(K1/2));
K3=h*f(T+(h/2),Y+(K2/2));
K4=h*f(T+h,Y+K3);
Y=Y+((K1+(2*(K2+K3))+K4)/6);
T=T+h;
%evaluation of RK4 at n=8192
end
A=Y;
%our true solution
for i=1:11
ERK2(i)=log10(abs(RK2(i)-A));
EMRK2(i)=log10(abs(MRK2(i)-A));
EEM(i)=log10(abs(EM(i)-A));
ERK4(i)=log10(abs(RK4(i)-A));
%calculating the log of error at T=2 for each n and each method
end
for i=1:11
LOGFE1(i)=log10(2*2^(i+1));
LOGFE2(i)=log10(2*2^(i+1));
LOGFE3(i)=log10(2^(i+1));
LOGFE4(i)=log10(4*2^(i+1));
%the log of the number of function evaluations
end
hold on
figure(1)
plot(LOGFE1,ERK2,'r')
plot(LOGFE2,EMRK2,'b')
plot(LOGFE3,EEM,'g')
plot(LOGFE4,ERK4,'y')
%plotting the logs of the errors with respect to log of function evalutions
hold off
%labeling all of our graphs
figure(1)
legend('RK2 error','My own RK2 error','EM error','RK4 error')
lgd=legend;
lgd.Title.String = 'Log of Errors';


format long
T=0;
Y=4;
r(1:13,1)=4;
%initial values
for j=3:15
n=2^(j);
%my number of steps
h=6/n;
%step size
for i=1:n
   K1=h*f(T,Y);
   K2=h*f(T+(h/2),Y+(K1/2));
   K3=h*f(T+(h/2),Y+(K2/2));
   K4=h*f(T+h,Y+K3);
   Y=Y+((K1+(2*(K2+K3))+K4)/6);
   T=T+h;
   r(j-2,i+1)=Y;
   %evaluation of RK4 using our fixed point function as the f inside our RK4
end
T=0;
Y=4;
end
hold on
for i=1:13
   %plotting each solution up to n=32768 steps
 p=[0:6/(2^(16-i)):6];
 plot(p,r(14-i,[1:((2^(16-i))+1)]),'color',[0 0 (i^3)/2197])
end
hold off
%labeling all of our graphs
figure(1)
legend('n=32768 steps')
%labeling the most accurate solution
lgd=legend;
lgd.Title.String = 'Solution';
function z=f(t,y)
%my function for root finding using fixed point
z=2;
for i=1:50
 z=(exp(-1-sin(z)))-(((sin(t + y))^2)*((1+(z^2))^(1/3)));
end
end