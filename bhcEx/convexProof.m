clear
syms mu k u q j J k0

OR=k0*sqrt(2);
OD=(q^j+1)*k0/sqrt(2);
OC=q^j*k0*sqrt(2);
OQ=(q^(J+1)+1)*k0/sqrt(2);
OG=(1+q^(-j))*q^(J+1)*k0/sqrt(2);
OS=q^(J+1)*k0*sqrt(2);

f = 0.5*(mu^2-k^2)-u*k^2;

intf1 = int(f,k,OC-mu,mu-OR);
fact = simplify(intf1/(2*mu-OR-OC));
a=simplify(diff(fact,mu,2)/2)
b=simplify(subs(diff(fact,mu,1),mu,0))
c=simplify(subs(fact,mu,0))
argmin=simplify(-b/2/a)
min1=simplify(subs(intf1,mu,OD))
max1=simplify(subs(intf1,mu,OC))

intf2 = int(f,k,mu-OC,mu-OR);
fact = simplify(intf2/(OC-OR));
a=simplify(diff(fact,mu,2)/2)
b=simplify(subs(diff(fact,mu,1),mu,0))
c=simplify(subs(fact,mu,0))
diff_intf2 = diff(intf2,mu);
argmax=simplify(-b/2/a)
min2=simplify(subs(intf2,mu,OC))
max2=simplify(subs(intf2,mu,OQ))

intf3 = int(f,k,mu-OC,OS-mu);
fact = simplify(intf3/(OS+OC-2*mu));
a=simplify(diff(fact,mu,2)/2)
b=simplify(subs(diff(fact,mu,1),mu,0))
c=simplify(subs(fact,mu,0))
argmin3=simplify(-b/2/a);
min3=simplify(subs(intf3,mu,OQ))
max3=simplify(subs(intf3,mu,OG))


