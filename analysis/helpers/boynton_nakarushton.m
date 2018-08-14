function R = boynton_nakarushton(C,a,p,q,sigma)

R = a*(C.^(p+q))./(C.^q+sigma^q);