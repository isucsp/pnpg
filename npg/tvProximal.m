function proximal=tvProximal(tv, project)
    % make an object to solve 0.5*||x-a||_2^2+u*TV(x)+I_C(x)
    % TV and C are specified via constructor see denoise for a and u
    proximal.penalty = @(z) tlv(z,tv);
    proximal.iterative=true;
    proximal.op=@denoise;

    if(~exist('project','var') || isempty(project))
        project=@(x)x;
    end

    prnt = 0;

    function [D,iter,pOut,fun_all]=denoise(Xobs,lambda,epsilon,MAXITER,pInit)
        % modified by renliang gu
        %
        %This function implements the FISTA method for TV-based denoising problems
        %
        % Based on the paper
        % Amir Beck and Marc Teboulle, "Fast Gradient-Based Algorithms for Constrained
        % Total Variation Image Denoising and Deblurring Problems"
        % -----------------------------------------------------------------------
        % Copyright (2008): Amir Beck and Marc Teboulle
        % 
        % FISTA is distributed under the terms of 
        % the GNU General Public License 2.0.
        % 
        % Permission to use, copy, modify, and distribute this software for
        % any purpose without fee is hereby granted, provided that this entire
        % notice is included in all copies of any software which is or includes
        % a copy or modification of this software and in all copies of the
        % supporting documentation for such software.
        % This software is being provided "as is", without any express or
        % implied warranty.  In particular, the authors do not make any
        % representation or warranty of any kind concerning the merchantability
        % of this software or its fitness for any particular purpose."
        % ----------------------------------------------------------------------
        % INPUT
        % Xobs ..............................an observed noisy image.
        % lambda ........................ parameter
        %  
        % OUTPUT
        % D ........................... The solution of the problem 
        %                                            min{||X-Xobs||^2+2*lambda*TV(X
        %                                            ) : l <= X_{i,j} <=u} 
        % iter .............................  Number of iterations required to get
        %                                            an optimal solution (up to a tolerance)
        % fun_all ......................   An array containing all the function
        %                                             values obtained during the
        %                                             iterations

        %Define the Projection onto the box

        if (~exist('MAXITER','var'))
            MAXITER=100;
        end
        if (~exist('epsilon','var'))
            epsilon=1e-4;
        end

        [m,n]=size(Xobs);
        D=zeros(m,n);
        if(exist('pInit','var') && ~isempty(pInit))
            P1=pInit{1}; P2=pInit{2};
        else
            [P1,P2]=Ltrans(D);
        end

        switch tv
            case 'iso'
                A=[P1;zeros(1,n)].^2+[P2,zeros(m,1)].^2;
                A=sqrt(max(A,1));
                P1=P1./A(1:m-1,:);
                P2=P2./A(:,1:n-1);
            case {'l1', '1d'}
                P1=max(min(P1,1),-1);
                P2=max(min(P2,1),-1);
            otherwise
                error('unknown type of total variation. should be iso or l1');
        end
        R1=P1; R2=P2;

        % debug=true;
        f=@(x,u) 0.5*norm(x-Xobs,'fro')^2+u*tlv(x,tv);

        tk=1;
        tkp1=1;
        count=0;
        i=0;

        fval=inf;
        fun_all=[];
        if(prnt)
            fprintf('***********************************\n');
            fprintf('*Solving with FGP/FISTA**\n');
            fprintf('***********************************\n');
            fprintf('#iteration  function-value  relative-difference\n');
            fprintf('---------------------------------------------------------------------------------------\n');
        end
        while((i<MAXITER)&&(count<5))
            fold=fval;
            %%%%%%%%%
            % updating the iteration counter
            i=i+1;
            %%%%%%%%%
            % Storing the old value of the current solution
            Dold=D;
            %%%%%%%%%%
            %Computing the gradient of the objective function
            Pold1=P1;
            Pold2=P2;
            tk=tkp1;
            C=Xobs-lambda*Lforward(R1,R2);
            D=project(C);
            [Q1,Q2]=Ltrans(D);

            %%%%%%%%%%
            % Taking a step towards minus of the gradient
            P1=R1+Q1/(8*lambda);
            P2=R2+Q2/(8*lambda);

            %%%%%%%%%%
            % Peforming the projection step
            switch tv
                case 'iso'
                    A=[P1;zeros(1,n)].^2+[P2,zeros(m,1)].^2;
                    A=sqrt(max(A,1));
                    P1=P1./A(1:m-1,:);
                    P2=P2./A(:,1:n-1);
                case {'l1','1d'}
                    P1=max(min(P1,1),-1);
                    P2=max(min(P2,1),-1);
                    %     P1=P1./(max(abs(P1),1));
                    %     P2=P2./(max(abs(P2),1));
                otherwise
                    error('unknown type of total variation. should be iso or l1');
            end

            %%%%%%%%%%
            %Updating R and t
            tkp1=(1+sqrt(1+4*tk^2))/2;

            R1=P1+(tk-1)*(P1-Pold1)/tkp1;
            R2=P2+(tk-1)*(P2-Pold2)/tkp1;

            re=norm(D-Dold,'fro')/norm(D,'fro');
            if (re<epsilon)
                count=count+1;
            else
                count=0;
            end

            fval=(-norm(C-D,'fro')^2+norm(C,'fro')^2)/2;
            if (fval>fold)
                tkp1=1;
            end

            if(nargout==4)
                fun_all=[fun_all;fval];
                if(prnt)
                    fprintf('%7d %10.10g %10.10g  %19g',i,fval,norm(D-Dold,'fro')/norm(D,'fro'),f(D,lambda));
                    if (fval>fold) fprintf('  *\n'); else fprintf('   \n'); end
                end
            end
        end
        iter=i;
        pOut={P1, P2};
    end

end

% If edit the following, update TV.[A,B]t\=
function x = Psi_v(p)
    [I,~]=size(p);
    p(I,:)=0;
    x=[p(1,:); p(2:I,:)-p(1:I-1,:)];
end
function p = Psi_vt(x)
    [I,J]=size(x);
    p=[x(1:I-1,:)-x(2:I,:);zeros(1,J)];
end
function x = Psi_h(q)
    [~,J]=size(q);
    q(:,J)=0;
    x=[q(:,1), q(:,2:J)-q(:,1:J-1)];
end
function q = Psi_ht(x)
    [I,J]=size(x);
    q=[x(:,1:J-1)-x(:,2:J), zeros(I,1)];
end

