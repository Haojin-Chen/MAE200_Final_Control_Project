function dX = PRiccat(E,A,C,P,Q2,Q1)
    dP_pre = A'*P*E + E'*P*A - E'*P*C'*inv(Q2)*C*P*E + Q1;
    dX = inv(E') * dP_pre * inv(E);
end
