for i=size(A,3):-1:2
        f1 = RHS(X(:, :, i), E(:, :, i), A(:, :, i), BC(:, :, i), R, Q);
        f2 = RHS(X(:, :, i)-f1/2h, E(:, :, i), A(:, :, i), BC(:, :, i), R, Q);
        f3 = RHS(X(:, :, i)-f2/2h, E(:, :, i), A(:, :, i), BC(:, :, i), R, Q);
        f4 = RHS(X(:, :, i)-f3h, E(:, :, i), A(:, :, i), BC(:, :, i), R, Q);
        X(:, :, i-1) = X(:, :, i) - h(f1/6+(f2+f3)/3+f4/6);
%     f1=-inv(transpose(E(:,:,i)))transpose(A(:,:,i))X(:,:,i)-X(:,:,i)A(:,:,i)inv(E(:,:,i))+X(:,:,i)BC(:,:,i)inv(R)transpose(BC(:,:,i))X(:,:,i)-inv(transpose(E(:,:,i)))QE(:,:,i);
%     f2=-inv(transpose(E(:,:,i)))transpose(A(:,:,i))(X(:,:,i)-hf1/2)-(X(:,:,i)-hf1/2)A(:,:,i)inv(E(:,:,i))+(X(:,:,i)-hf1/2)BC(:,:,i)inv(R)transpose(BC(:,:,i))(X(:,:,i)-hf1/2)-inv(transpose(E(:,:,i)))Qinv(E(:,:,i));
%     f3=-inv(transpose(E(:,:,i)))transpose(A(:,:,i))(X(:,:,i)-hf2/2)-(X(:,:,i)-hf2/2)A(:,:,i)inv(E(:,:,i))+(X(:,:,i)-hf2/2)BC(:,:,i)inv(R)transpose(BC(:,:,i))(X(:,:,i)-hf2/2)-inv(transpose(E(:,:,i)))Qinv(E(:,:,i));
%     f4=-inv(transpose(E(:,:,i)))transpose(A(:,:,i))(X(:,:,i)-hf3)-(X(:,:,i)-hf3)A(:,:,i)inv(E(:,:,i))+(X(:,:,i)-hf3)BC(:,:,i)inv(R)transpose(BC(:,:,i))(X(:,:,i)-hf3)-inv(transpose(E(:,:,i)))Qinv(E(:,:,i));
%     X(:,:,i-1)=X(:,:,i)-h(f1/6+(f2+f3)/3+f4/6);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dX=RHS(X, E,A,B,R,Q)
% keyboard
dX = - (E')^-1 *A' X - X * A * (E)^-1 + X * B * R^-1 * B' * X - (E')^-1 * Q * E^-1;
end


function R=RHS_time_march(x,u,s)
E=Compute_E(x,s); N=Compute_N(x,u,s); R=E\N;
end % function RHS
function N=Compute_N(x,u,s)
N=[x(4); x(5); x(6); -s.m1s.ell1sin(x(2))x(5)^2-s.m2s.ell2sin(x(3))x(6)^2+u;
 s.m19.8s.ell1sin(x(2)); s.m29.8s.ell2sin(x(3)) ];
end % function Compute_N

function R=RHS_NR(x,s);  g=9.8;
a42=s.m1s.ell1(x(8)sin(x(2))+x(5)^2cos(x(2))); a45=2s.m1s.ell1x(5)sin(x(2));
a43=s.m2s.ell2(x(9)sin(x(3))+x(6)^2cos(x(3))); a46=2s.m2s.ell2x(6)sin(x(3));
a52=s.m1s.ell1(gcos(x(2))-x(7)sin(x(2))); a63=s.m2s.ell2(gcos(x(3))-x(7)sin(x(3)));
A=[zeros(3) eye(3); 0 -a42 -a43 0 -a45 -a46; 0 a52 0 0 0 0; 0 0 a63 0 0 0];
R=A;
end % function RHS