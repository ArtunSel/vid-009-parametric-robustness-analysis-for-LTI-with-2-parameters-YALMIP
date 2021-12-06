clear all,close all,clc;
%%
% this is for SDPT3
current_dir=pwd;
cd('C:\Users\Artun\Desktop\BACK UP 2021-june 20\Documents 2020\MATLAB toolboxes for optimization\SDPT3-4.0');
run('Installmex.m')
run('startup.m')
cd(current_dir);

% this is for YALMIP
addpath(genpath('C:\Users\Artun\Desktop\BACK UP 2021-june 20\Documents 2020\MATLAB toolboxes for optimization\YALMIP-master'))
yalmip('clear');
%% yalmip param robust ness anlaysis
clear all,close all,clc;yalmip('clear');

p1=3; % takes values b/w [0,9]
p2=7; % takes values b/w [2,10]

A=[-p1    ,1,-1;
    2-2*p2,2,-1;
    3     ,1,-p1*p2];
% Que.: does there exist a P=P'>0 s.t. "<P*A>_s < 0"
eps1=0.01;
P=sdpvar(3,3,'symmetric');
F=[];
F=[F;P>=eps1*eye(3)];
F=[F;A'*P+P*A<=-eps1*eye(3)];
F=[F;-10<=vec(P)<=10];

ops = sdpsettings('solver','sdpt3');
sol = optimize(F,[],ops);
sol.info

P0=value(P)
min(eig(P0))
eig(A'*P0+P0*A)
%% TASK-1
tic %-----------------------------------------------
discretization_num=30;
p1_vec=linspace(0,9 ,discretization_num);
p2_vec=linspace(2,10,discretization_num);
blue_dot_locations=[];

P=sdpvar(3,3,'symmetric');
for ii=1:1:length(p1_vec)
    for jj=1:1:length(p2_vec)
        temp_p1=p1_vec(ii);
        temp_p2=p2_vec(jj);
        %
        A=[-temp_p1    ,1,-1;
            2-2*temp_p2,2,-1;
            3     ,1,-temp_p1*temp_p2];
        % YALMIP feasiblity problem
        eps1=0.01;
        F=[];
        F=[F;P>=eps1*eye(3)];
        F=[F;A'*P+P*A<=-eps1*eye(3)];
        
        ops = sdpsettings('solver','sdpt3');
        sol = optimize(F,[],ops);
%         sol.info
        if(sol.problem==1) % that'd indicate the "infeasibility"
%             plot(temp_p1,temp_p2,'r*'); hold on;
        else
%             plot(temp_p1,temp_p2,'b*'); hold on;
            blue_dot_locations=[blue_dot_locations;[temp_p1,temp_p2]];
        end
    end
end
toc %-----------------------------------------------
fig1=figure(1);fig1.Color=[1,1,1];
plot(blue_dot_locations(:,1),blue_dot_locations(:,2),'b*');
xlabel('p1');ylabel('p2');
%% QUADRATIC STABILITY TASK-2 prove Q.S. for a given param. intervals

A0 =[0,1,-1;2,2,-1;3,1,1];
Ap1=zeros(3); Ap1(1,1)=-1;
Ap2=zeros(3); Ap2(2,1)=-2;
Ap3=zeros(3); Ap3(3,3)=-1;

p1_min=3-0.2;   p1_max=3+0.2;
p2_min=7-0.5;   p2_max=7+0.5;
p3_min=p1_min*p2_min;   p3_max=p1_max*p2_max;

p1_val=[p1_min,p1_max];
p2_val=[p2_min,p2_max];
p3_val=[p3_min,p3_max];
% construct the A_i matrices i=1...8 [2x2x2]
A_1=A0+p1_val(1+ 0 )*Ap1+p2_val(1+ 0 )*Ap2+p3_val(1+ 0 )*Ap3;     % 0,0,0
A_2=A0+p1_val(1+ 1 )*Ap1+p2_val(1+ 0 )*Ap2+p3_val(1+ 0 )*Ap3;     % 1,0,0
A_3=A0+p1_val(1+ 0 )*Ap1+p2_val(1+ 1 )*Ap2+p3_val(1+ 0 )*Ap3;     % 0,1,0
A_4=A0+p1_val(1+ 1 )*Ap1+p2_val(1+ 1 )*Ap2+p3_val(1+ 0 )*Ap3;     % 1,1,0

A_5=A0+p1_val(1+ 0 )*Ap1+p2_val(1+ 0 )*Ap2+p3_val(1+ 1 )*Ap3;     % 0,0,1
A_6=A0+p1_val(1+ 1 )*Ap1+p2_val(1+ 0 )*Ap2+p3_val(1+ 1 )*Ap3;     % 1,0,1
A_7=A0+p1_val(1+ 0 )*Ap1+p2_val(1+ 1 )*Ap2+p3_val(1+ 1 )*Ap3;     % 0,1,1
A_8=A0+p1_val(1+ 1 )*Ap1+p2_val(1+ 1 )*Ap2+p3_val(1+ 1 )*Ap3;     % 1,1,1


eps1=0.01;
P=sdpvar(3,3,'symmetric');
F=[];
F=[F;P>=eps1*eye(3)];
F=[F;A_1'*P+P*A_1<=-eps1*eye(3)];
F=[F;A_2'*P+P*A_2<=-eps1*eye(3)];
F=[F;A_3'*P+P*A_3<=-eps1*eye(3)];
F=[F;A_4'*P+P*A_4<=-eps1*eye(3)];

F=[F;A_5'*P+P*A_5<=-eps1*eye(3)];
F=[F;A_6'*P+P*A_6<=-eps1*eye(3)];
F=[F;A_7'*P+P*A_7<=-eps1*eye(3)];
F=[F;A_8'*P+P*A_8<=-eps1*eye(3)];
F=[F;-50<=vec(P)<=50];


ops = sdpsettings('solver','sdpt3');
sol = optimize(F,[],ops);
sol.info
%%
P0=value(P)
min(eig(P0))
max(eig(A_1'*P0+P0*A_1))
max(eig(A_2'*P0+P0*A_2))
max(eig(A_3'*P0+P0*A_3))
max(eig(A_4'*P0+P0*A_4))

max(eig(A_5'*P0+P0*A_5))
max(eig(A_6'*P0+P0*A_6))
max(eig(A_7'*P0+P0*A_7))
max(eig(A_8'*P0+P0*A_8))






















%





