
%% (A) Exp.
u;
y;

plot(t,y(1:2:end))
%% (B) LTI

Y_lti = H*U
plot(t,Y_lti(1:2:end))

%% (C) nonlin. model

% Full input, org. out
[Ad_con,Bd_con,Cd_con,Dd_con] = systemMatriciesSS_dis(M,K,C,dof,in_dof_con,out_dof,opt.out_type,dt);
[H_FI] = ToeplitzMatrix(N,ms,r_con,Ad_con,Bd_con,Cd_con,Dd_con);

Gamma = pinv(H_FI)*(Y - H*U);

Y_est = H*U + H_FI*Gamma;

y_est = reshape(Y_est, ms, N);  % decollapse dof columns

plot(t,Y_est(1:2:end))


%% (D) include extended output

% H_ex: full output, org. input


% full input, full output
H_FIFO = H_con;

Y_estex = H_ex*U + H_FIFO*Gamma






%%
figure
tiledlayout('flow')
for i = 1:4
    nexttile
    hold on
    plot(t,Y_acc(i:4:end))
    plot(t,Y_estex(i:4:end))
    legend('actual','est.')
    title(sprintf('DOF: %d', i));

end




