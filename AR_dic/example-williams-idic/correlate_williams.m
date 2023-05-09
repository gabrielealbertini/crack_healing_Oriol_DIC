function B = correlate_williams(B)
maxit = B.maxit;
convcrit = B.convcrit;
init  = B.init;
options = B.options;
Mesh = B.Mesh;
Ip1 = B.Ip1;
Ip2 = B.Ip2;
Im1 = B.Im1;
scal = B.scal;
trustregion = B.trustregion;
dampen = B.dampen;
alpha = B.alpha;
cstate = B.cstate;
a = B.a;
crackpath = B.crackpath;
f = B.f;
g = B.g;


fprintf('%3s, %10s, %10s, %10s, %10s, %10s, %10s \n','it','a','da','K1','K2','TS','r');
for it = 1:maxit
    
    % Specifying the Sensitivity matrix for the current parameters.
    options.S = Williams_S(Mesh);
    
    % Integrated DIC
    % -----------------
    
    % Correlate with the Integrated DIC sensitivity fields
    cor = correlate_2D(f,g,Mesh,init,options);
    
    init = cor.U;
    
    % Stress Intensitiy factors:
    K1 = cor.U(Ip1(1))*scal(Ip1(1));
    K2 = cor.U(Ip1(2))*scal(Ip1(2));
    
    % T-Stress
    TS = 4*cor.U(Ip2(1))*scal(Ip2(1))/sqrt(2*pi);
    
    % Estimated crack tip shift
    if abs(K1) > abs(K2)
        % Mode I dominant
        da = -2*Mesh.Rnorm*cor.U(Im1(1))/cor.U(Ip1(1));
    else
        % Mode II dominant
        da = -2*Mesh.Rnorm*cor.U(Im1(2))/cor.U(Ip1(2));
    end
    
    % Trick 1: if the crack tip shift is too large, then don't trust it
    if trustregion > 0 % set to true to enable
        if da > trustregion
            da = trustregion;
        elseif da < -trustregion
            da = -trustregion;
        end
    end
    
    % Trick 2: set a dampening ratio for the NR scheme (0 < dampen < 2)
    if dampen > 0 % set to true to enable
        da = dampen * da;
    end
    
    % Trick 3: apply an exponential moving average to stabilize over the
    % iterations.
    if alpha > 0 % set to true to enable
        % 0 < alpha <= 1, the smaller alpha the longer the moving average
        if it == 1
            % initialze the average
            avg_a = a;
        end
        % update the moving average
        avg_a = (1-alpha)*(a+da) + alpha*avg_a;
        if it > 5
            % only use the avg after a few iterations
            da = avg_a - a;
        end
    end
    
    % move the crack tip Mesh to the new position
    Mesh = crack_movemesh(crackpath,Mesh,a,da);
    
    % update the crack length
    a = a + da;
    
    % update the crack position
    crackpos = crack_position(crackpath,a);
    crackang = crack_angle(crackpath,a);
    
    % update the stored parameters (for the Williams_S function)
    Mesh.crackpos = crackpos;
    Mesh.crackang = crackang;
    
    % print a status line (red if not converged properly)
    fid = (cor.s(1) ~= 0) + 1;
    fprintf(fid,'%3d, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e %s \n',it,a,da,K1,K2,TS,cor.r,cstate{cor.s(1)+1});
    
    % Convergence testing
    if abs(da)<convcrit
        break
    end
    
end

% add results to the input/output structure
B.init = cor.U;
B.cor = cor;
B.Mesh = Mesh;
B.a = a;
B.da = da;
B.it = it;
B.K1 = K1;
B.K2 = K2;
B.TS = TS;




