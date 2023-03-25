function dX = KRiccat(E,A,B,X,R,Q)
    dX_pre = A'*X*E + E'*X*A -E'*X*B*inv(R)*B'*X*E + Q
    dX = -inv(E') * dX_pre * inv(E)
end
