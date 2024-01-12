parten = @(x) sign(x).*(abs(x)-floor(abs(x)));

# Función de la familia de Arnold y cómputo de la melodía.

alpha = e;
eepsilon = 0.5;

phi = @(x) x+alpha+(eepsilon/(2*pi))*sin(2*pi*x);
P0 = [0];
for k = 1:15 P0 = [P0 parten(phi(P0(end)))]; end

# Cómputo de la primera aproximación al cambio de coordenadas y primera variación de la melodía.

P1 = [0];
H = @(x) x - (eepsilon/(2*pi))*cos(2*pi*x-alpha*pi)/(2*sin(alpha*pi));
for k = 1:15
 u = phi(H(P1(end)));
 F = @(x) H(x)-u;
 v = mod(fsolve(F,phi(P1(end))),1);
 P1 = [P1 v];
end

# Cómputo de la segunda aproximación al cambio de coordenadas y segunda variación de la melodía.

P2 = [0];
H = @(x) x - (eepsilon/(2*pi))*cos(2*pi*x-alpha*pi)/(2*sin(alpha*pi))+(eepsilon^2/(4*pi^2))*sin(4*pi*x-alpha*pi)/(4*sin(alpha*pi)*sin(2*alpha*pi));
for k = 1:15
 u = phi(H(P2(end)));
 F = @(x) H(x)-u;
 v = mod(fsolve(F,phi(P2(end))),1);
 P2 = [P2 v];
end

# Función de cambio de coordenadas en términos de los coeficientes de Fourier de la perturbación eta.

function resp = hache(coeffp,coeffn,alpha,x)
 m = rows(coeffp);
 n = columns(coeffp);
 aux = x;
 for l = 1:m
  upos=coeffp(l,:)./(exp(2*pi*i*[1:n]*alpha)-1);
  uneg=coeffn(l,:)./(exp(2*pi*i*-[1:n]*alpha)-1);
  resp = [];
  for k = 1:length(aux)
   resp = [resp sum(upos.*exp(2*pi*i*[1:n]*aux(k)))+sum(uneg.*exp(2*pi*i*-[1:n]*aux(k)))];
  end
  aux = aux + resp;
 end
 resp = aux;
end

function resp = evdif(mcoeffp,mcoeffn,phi,alpha,x)
 q = rows(mcoeffp);
 resp = x;
  u = phi(hache(mcoeffp,mcoeffn,alpha,resp));
  F = @(t) hache(mcoeffp,mcoeffn,alpha,t)-u;
  resp = fsolve(F,resp+alpha);
end

function P = kam_generator(alpha,epsilon,n,M)
 # Generates an array of M rows and N colums where the first row
 # corresponds to a melody generated with Arnold's diffeomorphism
 # with parameters alpha and epsilon, and M-1 changes of coordinates
 # and its corresponding variations of the original melody.

 # Function of the Arnold's family.

 phi = @(x) x+alpha+(epsilon/(2*pi))*sin(2*pi*x);

 NF = 10; # Number of Fourier coefficients.
 NP = 64; # Number of sampling points.
 
 t = linspace(0,1,NP+1); # Sampled points in the unit interval.
 
 coeffp = zeros(1,NF);
 coeffn = zeros(1,NF);

 P = zeros(M,n);

 for j = 1:M
  
  # Computation of the melody.
  
  for k = 2:n
   P(j,k) = mod(evdif(coeffp,coeffn,phi,alpha,P(j,k-1)),1);
  end
  
  ccoeffp = zeros(1,NF);
  ccoeffn = zeros(1,NF);
  res = [];
  
  #  rho = evdif(coeffp,coeffn,phi,alpha,evdif(coeffp,coeffn,phi,alpha,evdif(coeffp,coeffn,phi,alpha,0)))/3
 
  # Computation of eta in the sampled points.
 
  for q = t(1:end)
   v = evdif(coeffp,coeffn,phi,alpha,q)-q-alpha;
   res = [res v];
  end
 
  # Computation of the Fourier coefficients of eta.
  
  for s = 1:NF
   ccoeffp(s)=sum(res(2:end).*exp(-2*pi*i*s*t(2:end)))*(t(2)-t(1));
  end
  for s = 1:NF
   ccoeffn(s)=sum(res(2:end).*exp(2*pi*i*s*t(2:end)))*(t(2)-t(1));
  end
  coeffp = [coeffp; ccoeffp];
  coeffn = [coeffn; ccoeffn];
 end
end