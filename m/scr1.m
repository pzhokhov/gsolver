
load const_SI light_v vacuum_permittivity e e_mass

n2 = 3.45e-20;
omega0 = 2*pi*light_v/800e-9;
roc = 1.4533^2*vacuum_permittivity*e_mass*omega0^2/e/e

tau = [0.4 1 3];
E = (2:2:30)';

mI = zeros(size(E*tau));
mro = mI;
mdn = mI; 

for k = 1:size(tau,2);
for i=1:size(E,1);
     if k==1
        s = sprintf('tb%dnJ_av04.dump.parsed',E(i));
     else 
        s = sprintf('tb%dnJ_av%d.dump.parsed',E(i),tau(k));
     end;
     n = dir(s);
     if isempty(n) 
         warning(sprintf('File %s not found!',s));
         continue;
     end;
     [ro, I, f] = load_parseddump(s);
     mro(i,k) = max(ro(:));
     mI(i,k)  = max(I(:));
end;
mdn(:,k) = n2*mI(:,k) + 1.4533*(real(sqrt(1-mro(:,k)./roc/(1+j/(omega0*tau(k))))-1)); 
end;



