function [D,C,G,B,K,dV,dVl,Al,Bl,H,LfH,dLfH]=dynamics_three_link(x,a)
% DYNAMICS_THREE_LINK    Model of three-link biped walker model.
%    [D,C,G,B,K,dV,dVl,Al,Bl,H,LfH,DLFH] = DYNAMICS_THREE_LINK(X,
%      A) is the three-link
%    biped walking model. (x is of dimension 6)

% Eric Westervelt
% 23-Aug-2022 17:42:54

[r,m,Mh,Mt,L,g]=model_params_three_link;

[th3d,th1d,alpha,epsilon]=control_params_three_link;

th1=x(1); th2=x(2); th3=x(3);
dth1=x(4); dth2=x(5); dth3=x(6);

% D matrix
D=zeros(3);
D(1,1)=Mh*r^2 + Mt*r^2 + (5*m*r^2)/4;
D(1,2)=-(m*r^2*cos(th1 - th2))/2;
D(1,3)=L*Mt*r*cos(th1 - th3);
D(2,1)=-(m*r^2*cos(th1 - th2))/2;
D(2,2)=(m*r^2)/4;
D(3,1)=L*Mt*r*cos(th1 - th3);
D(3,3)=L^2*Mt;

% C matrix
C=zeros(3);
C(1,2)=-(dth2*m*r^2*sin(th1 - th2))/2;
C(1,3)=L*Mt*dth3*r*sin(th1 - th3);
C(2,1)=(dth1*m*r^2*sin(th1 - th2))/2;
C(3,1)=-L*Mt*dth1*r*sin(th1 - th3);

% G matrix
G=zeros(3,1);
G(1)=- Mh*g*r*sin(th1) - Mt*g*r*sin(th1) - (3*g*m*r*sin(th1))/2;
G(2)=(g*m*r*sin(th2))/2;
G(3)=-L*Mt*g*sin(th3);

% B matrix
B=zeros(3,2);
B(1,1)=-1;
B(2,2)=-1;
B(3,1)=1;
B(3,2)=1;

% K matrix
K=zeros(2,4);
K(1,1)=156;
K(1,3)=25;
K(2,2)=110;
K(2,4)=21;

% dV matrix
dV=zeros(1,6);
dV(1,1)=(184*th1)/77 - 10*dth2 - (242673938205*dth3)/38685626227668133590597632 - 10*dth1 + (184*th2)/77 + (60052912102237*th3)/79228162514264337593543950336 - (60052912102237*th3d)/158456325028528675187087900672 - (60052912102237*conj(th3d))/158456325028528675187087900672;
dV(1,2)=(184*th1)/77 - 10*dth2 - (242673938205*dth3)/38685626227668133590597632 - 10*dth1 + (184*th2)/77 + (60052912102237*th3)/79228162514264337593543950336 - (60052912102237*th3d)/158456325028528675187087900672 - (60052912102237*conj(th3d))/158456325028528675187087900672;
dV(1,3)=(60052912102237*th1)/79228162514264337593543950336 - (25836402041767*dth2)/9903520314283042199192993792 - 10*dth3 - (25836402041767*dth1)/9903520314283042199192993792 + (60052912102237*th2)/79228162514264337593543950336 + (391*th3)/195 - (391*th3d)/390 - (391*conj(th3d))/390;
dV(1,4)=(370*dth1)/7 + (370*dth2)/7 + (8848846829079*dth3)/309485009821345068724781056 - 10*th1 - 10*th2 - (25836402041767*th3)/9903520314283042199192993792 + (25836402041767*th3d)/19807040628566084398385987584 + (25836402041767*conj(th3d))/19807040628566084398385987584;
dV(1,5)=(370*dth1)/7 + (370*dth2)/7 + (8848846829079*dth3)/309485009821345068724781056 - 10*th1 - 10*th2 - (25836402041767*th3)/9903520314283042199192993792 + (25836402041767*th3d)/19807040628566084398385987584 + (25836402041767*conj(th3d))/19807040628566084398385987584;
dV(1,6)=(8848846829079*dth1)/309485009821345068724781056 + (8848846829079*dth2)/309485009821345068724781056 + (314*dth3)/5 - (242673938205*th1)/38685626227668133590597632 - (242673938205*th2)/38685626227668133590597632 - 10*th3 + 5*th3d + 5*conj(th3d);

% dVl matrix
dVl=zeros(1,4);
dVl(1,1)=(60052912102237*th1)/79228162514264337593543950336 - (25836402041767*dth2)/9903520314283042199192993792 - 10*dth3 - (25836402041767*dth1)/9903520314283042199192993792 + (60052912102237*th2)/79228162514264337593543950336 + (391*th3)/195 - (391*conj(th3d))/195;
dVl(1,2)=(184*th1)/77 - 10*dth2 - (242673938205*dth3)/38685626227668133590597632 - 10*dth1 + (184*th2)/77 + (60052912102237*th3)/79228162514264337593543950336 - (60052912102237*conj(th3d))/79228162514264337593543950336;
dVl(1,3)=(8848846829079*dth1)/309485009821345068724781056 + (8848846829079*dth2)/309485009821345068724781056 + (314*dth3)/5 - (242673938205*th1)/38685626227668133590597632 - (242673938205*th2)/38685626227668133590597632 - 10*th3 + 10*conj(th3d);
dVl(1,4)=(370*dth1)/7 + (370*dth2)/7 + (8848846829079*dth3)/309485009821345068724781056 - 10*th1 - 10*th2 - (25836402041767*th3)/9903520314283042199192993792 + (25836402041767*conj(th3d))/9903520314283042199192993792;

% Al matrix
Al=zeros(4,4);
Al(1,3)=1;
Al(2,4)=1;

% Bl matrix
Bl=zeros(4,2);
Bl(3,1)=1;
Bl(4,2)=1;
a01=a(1); a11=a(2); a21=a(3); a31=a(4);
a02=a(5); a12=a(6); a22=a(7); a32=a(8);

% Ha matrix
H=zeros(2,1);
H(1,1)=th3 - a01 - a11*th1 - a21*th1^2 - a31*th1^3;
H(2,1)=th1 + th2 - (th1 + th1d)*(th1 - th1d)*(a02 + a12*th1 + a22*th1^2 + a32*th1^3);

% LfHa matrix
LfH=zeros(2,1);
LfH(1,1)=dth3 - dth1*(a11 + 2*a21*th1 + 3*a31*th1^2);
LfH(2,1)=dth2 - dth1*((th1 - th1d)*(a02 + a12*th1 + a22*th1^2 + a32*th1^3) + (th1 + th1d)*(a02 + a12*th1 + a22*th1^2 + a32*th1^3) + (th1 + th1d)*(th1 - th1d)*(a12 + 2*a22*th1 + 3*a32*th1^2) - 1);

% dLfHa matrix
dLfH=zeros(2,6);
dLfH(1,1)=-dth1*(2*a21 + 6*a31*th1);
dLfH(1,4)=- a11 - 2*a21*th1 - 3*a31*th1^2;
dLfH(1,6)=1;
dLfH(2,1)=-dth1*(2*a02 + 2*(th1 + th1d)*(a12 + 2*a22*th1 + 3*a32*th1^2) + 2*a12*th1 + 2*(th1 - th1d)*(a12 + 2*a22*th1 + 3*a32*th1^2) + 2*a22*th1^2 + 2*a32*th1^3 + (th1 + th1d)*(2*a22 + 6*a32*th1)*(th1 - th1d));
dLfH(2,4)=1 - (th1 + th1d)*(a02 + a12*th1 + a22*th1^2 + a32*th1^3) - (th1 + th1d)*(th1 - th1d)*(a12 + 2*a22*th1 + 3*a32*th1^2) - (th1 - th1d)*(a02 + a12*th1 + a22*th1^2 + a32*th1^3);
dLfH(2,5)=1;