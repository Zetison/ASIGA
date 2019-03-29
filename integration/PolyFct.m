function [W,G,U] = PolyFct(p)
 %
 %      *** Factoring a multiple-root polynomial ***
 %
 %     A given polynomial is factored into lower-degree
 %   distict-root polynomials with natual-order-integer
 %   powers. Roots with multiplicites may then be found
 %   with easy.
 %
 %     The more multiplicities the polynomial roots possess,
 %   the more efficient this routine will be.  This is
 %   contrary to the general public statement.
 %
 %     If all coefficients of given polynomial are integers,
 %   then the derived routine involves only very simple
 %   arithmetric operations, such as pure integer addition
 %   and multiplication. There are no floating point oprations
 %   in process. The round-off errors may thus be eliminated,
 %   iff there is no numerical digital overflow.
 %
 %     The crucial concern is polynomial GCD computation.
 %   The classical Euclidean GCD algorithm is applied here 
 %   for its simplicity in mathematics, however, it is not 
 %   stable numerically.
 %
 %     The routine would be very useful, if a suitable GCD
 %   computation algorithm could be found with application
 %   of any higher mathematical techniques.
 %
 %   Reference:
 %     F C Chang,"Factorization of Multiple-Root Polynomials" 
 %     via  fcchang007@yahoo.com
 %                                              03/18/2008
 %
 %  
 %   Description of MATLAB M-file: PolyFct.m
 %
 %    [W,G,U] = PolyFct(p)
 %       Main routine: from given polynomial p
 %       to get factored polynomials W, and
 %       related sub-polynomials G and U.
 %    [qc,rc,q,r] = PolyDiv(b,a)
 %       Polynomial division of b and a
 %       to get quotient qc and reminder rc.
 %       Note:  b = a*q + r, qc = c*q, rc = c*r
 %     g = PolyGcd(b,a)
 %       Polynomial GCD of given b and a
 %       to get g by repeating polynomial divisions .
 %    [q,v] = recast(p)
 %       Recast vector p by dividing its own gcd=v
 %       to get q
 %
 %    ***  ONLY PURE INTEGER ARITHMETIC OPERATIONS  ***
 %    ***  NO FLOWTING POINT CALCULATIONS INVOLVED  ***
 %
 %     Amazingly, this simple routine (less than 70 lines)
 %   gives the exact results for the test polynomials of 
 %   fairly high degrees, such as
 %
 %       p(x) = (x + 1)^50
 %       p(x) = (x^4 - 1)^25
 %       p(x) = (x - 123456789)^4
 %       p(x) = (123x + 456)^4
 %       p(x) = (x^4 -2x^3 +3x^2 -4x +5)^12
 %
 %   EXAMPLE:
 %
 %   For given test polynomial:
 %
 %    p(x) =   x^19 -3x^18 -8x^17 -4x^16 +76x^15 +284x^14
 %          -536x^13 -808x^12 -2474x^11 +7486x^10 +5896x^9
 %          -3872x^8 -2728x^7 -15812x^6 +77849x^5 +2184x^4
 %          -82319x^3 +24045x^2 +28800x -13500
 %
 %   we shall get factorization of original polynomial:
 %
 %    p(x) = (x +3)^1 (x^2 -5x +6)^2 (x^3 +3x^2 +7x +5)^3
 %           (1)^4 (x -1)^5
 %
 %   All roots with multiplicities may thus be determined.
 %
 %  >>
 %  >> % Create a test polynomial of degree 19:
 %  >> p=poly([ 1 1 1 1 1 -1 -1 -1 -1+2i -1+2i -1+2i  ...
 %                           -1-2i -1-2i -1-2i 2 2 3 3 -3])
 %     p =
 %           1         -3         -8         -4         76 
 %         284       -536       -808      -2474       7486
 %        5896      -3872     -27284     -15812      77848
 %        2184     -82319      24045      28800     -13500
 %
 %  >> % Get all factors of lower-degree distinct polynomials
 %  >> W = PolyFct(p); celldisp(W)
 %
 %     W{1} =
 %             1     3
 %     W{2} =
 %             1    -5     6
 %     W{3} =
 %             1     3     7     5
 %     W{4} =
 %             1
 %     W{5} =
 %             1    -1
 %
 %  >> % End of Example
 %  >>
 %

 
%function [W,G,U] = PolyFct(p)
 %    *** Factorization of a multiple-root polynomial ***
 %   By  F C Chang     3/18/2008
 %
      g2 = p;
  for k = 1:length(p);
      g1 = g2;
      g2 = PolyGcd(g1,polyder(g1));
      g3 = PolyGcd(g2,polyder(g2));
      u1 = PolyDiv(g1,g2);
      u2 = PolyDiv(g2,g3);
      w1 = PolyDiv(u1,u2);
      G{k} = g1;  U{k} = u1;  W{k} = w1;
    if length(g2) == 1;    break;   end; 
  end;
      G{k+1} = g2;  G{k+2} = 1;   U{k+1} = u2;
    % celldisp(G);  celldisp(U);  celldisp(W);

function g = PolyGcd(b,a)
 %   GCD of a pair of polynomials b and a
 %  
       p2 = b;   p3 = a;
   while(1);
       p1 = recast(p2);
       p2 = recast(p3);
      [qc,rc] = PolyDiv(p1,p2);
       p3 = recast(rc);
    if norm(p3,inf) < 1.e-6;  break;  end    
   end;
       g = p2;

function [qc,rc,q,r] = PolyDiv(b,a)
 %   Polynomial division by convolution matrix. 
 %
        n = length(b);   m = length(a);
       [b,vb] = recast(b);  [a,va] = recast(a);
  if m > n;  q = 0;   r = vb*b;  qc = 0;  rc = b;  return; end;
  if m == 1; q = (vb/va)*b; r = 0; qc = b; rc = 0; return; end;
        qk(1) = a(1)^(n-m+1);
  for k = 2:n+1;
     if k < n-m+3;  qh(k-1) = qk(k-1)/a(1);
        qk(k) = [b(k-1),-a(min(k-1,m):-1:2)]...
               *[qh(1),qh(max(2,k+1-m):k-1)]';
     else
        qk(k) = [b(k-1),-a(min(m,k-1):-1:k-1+m-n)]...
               *[qk(1),qk(max(k+1-m,2):n-m+2)]';
     end;       
  end;
        c = qk(1);
       [qc,vq] = recast(qk(2:n-m+2)); 
        rc = qk(n-m+3:n+1);
        nz = find(abs(rc)>0);
     if isempty(nz);  rc = 0; vr = 1;  else
       [rc,vr] = recast(rc(nz(1):end));
     end;
        q = (vq*vb/va/c)*qc;
        r = (vr*vb/c)*rc;
         
function [q,v] = recast(p)
 %   Recast a vector by dividing it with its own GCD.
 %
        v = norm(p,inf);
    if abs(p) == 0;  q = abs(p); v = 1;  return; end;
    if length(p) > 1;  
       for j = 1:length(p),  v = gcd(v,p(j));  end;
    end;
        q = p/v;
    if p(1) ~= 0;  q=q*sign(p(1));  end;



