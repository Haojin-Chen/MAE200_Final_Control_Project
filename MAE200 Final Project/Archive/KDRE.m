function X_t = KDRE(X_T, T, h)

for i = 1:T/h
    dX = KRiccati(E,A,B,X,R,Q)

    k1 = h*dX;
    k2 = h*KRiccati(E,A,B,X,R,Q)




