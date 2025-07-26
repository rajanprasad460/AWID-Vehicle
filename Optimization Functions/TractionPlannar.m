function out = TractionPlannar(nd_approach,A,B,lb,ub,x0)
%================ [PRIORITY BASED FORCE DISTRIBUTION] ==============
switch nd_approach
    case 1
        x = PriorityBasedDistribution(A,B,lb,ub);

        %============== Optimization starts here  =============================

        %========== [FMINCON FUNCTION ]  ==================================
    case 2
        %---------------  [Fmincon option OBJECTIVE FUNCTION]  ----------
        fun= @(x)x*x';
        x = fmincon(fun,x0,[],[],A,B,lb,ub);
        x = x';
        %======================================================================
        % =========== [QP PROBLEM FORMULATION BASED OPTIMIZATION ]============

    case 3
        gamma=0.001;
        Wv = diag(ones(1,length(B)),0);
        %         Wv = [1 0 0;0 1 0;0 0 10];
        Wu=diag(ones(1,length(x0)),0);
        H=A'*Wv*A+gamma*Wu;    f=(-B'*Wv*A)';
        A_bd=[diag(ones(1,length(x0)),0);diag(-1*ones(1,length(x0)),0)];
        bd=[ub;-lb];
        % ----------- Select the method ----------------------------------
        %     x = quadprog(H,f,A_bd,bd); % Matlab QUADPROG function
        %     x= MSTHESIS_OPT_BARRIER_FULL(H,f,A_bd,bd,x0); % Barrier Method Original
        x= MSTHESIS_OPT_BARRIER(H,f,A_bd,bd,x0); % Customized Barrier Method
        %     x= MSTHESIS_OPT_IP(H,f,A_bd,bd,x0); % Customized Interior Point Method
        %========================================================================
    case 4
        gamma=0.001;
        Wv=diag(ones(1,length(B)),0);
        %         Wv = [1 0 0;0 1 0;0 0 10];
        Wu=diag(ones(1,length(x0)),0);
        H=A'*Wv*A+gamma*Wu;    f=(-B'*Wv*A)';
        A_bd=[diag(ones(1,length(x0)),0);diag(-1*ones(1,length(x0)),0)];
        bd=[ub;-lb];
        % ----------- Select the method ----------------------------------
        x = quadprog(H,f,A_bd,bd); % Matlab QUADPROG function
        %             x= quadprog(H,f,[],[],eye(4),lb,ub,x0);

        %========================================================================
        % ===============  [LSQLIN FUNCTION ] =================================
    case 5
        %             C_ls=diag(0.75*ones(1,4));
        %             d_ls=0.01*ones(4,1);
        A_bd=[diag(ones(1,4),0);diag(-1*ones(1,4),0)];
        bd=[ub;-lb];
        %     x = lsqlin(C_ls,d_ls,[],[],A,B,lb,ub,x0);
        x = lsqlin(A,B,A_bd,bd);
end
%%%%%%------------------------------------------------------------------
%========================================================================
out = x;
end