function A = generateLaplaceOnSquare(m, bdrCond)
    if strcmpi(bdrCond, 'dirichlet')
        I = speye(m-1);
        e = ones(m-1,1);
        T = spdiags([e -4*e e],[-1 0 1],m-1,m-1);
        S = spdiags([e e],[-1 1],m-1,m-1);
        A = (kron(I,T) + kron(S,I));
    else
        B = [ones(m*m,1),ones(m*m,1),-4*ones(m*m,1),ones(m*m,1),ones(m*m,1)];
        d = [-m,-1,0,1,m];
        A = spdiags(B,d,m^2,m^2);
        %adapt corner nodes
        A(1,2) = 2;
        A(1,m+1) = 2;
        A(m*m,end-1) = 2;
        A(m*m,end-m) = 2;

        %adapt remaining boundary nodes
        for i = 2 : m
            A(i,i+m) = 2;
            row = m^2-i+1;
            col = row-m;
            A(row,col) = 2;
        end

        %enforcing auxilliary Dirichlet condition in corner point,
        %in order to have a unique solution (equivalent to fixing
        %the mean)
        if strcmpi(bdrCond, 'neumann')
            A(1,:) = 0;
            A(1,1) = 1;
        end
    end 
    A = -A;
end