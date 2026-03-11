within EMBSlib.Internal;
package MatrixFunctions
 function getTaylorFunction
  input Integer nr;
  input Integer nq;
  input Integer nc;
  input Real[nr,nc] M0;
  input Real[nr,nq,nc] M1;
  input Real[nq] q;
  output Real[nr,nc] M;
  //annotation(derivative=getTaylorFunctionDerivative, derivative(zeroDerivative=M0)=getTaylorFunctionDerivative_onlyM1);
  protected
   Real[nr,nc] aux;
 algorithm
   M := zeros(nr,nc);
   for i in 1:nq loop
     aux := M1[:,i,:]*q[i];
     M := M+aux;
   end for;
   M := M0 +M;
   annotation(smoothOrder=1);
 end getTaylorFunction;

 function getGrMatrix
  input Integer nr;
  input Integer nq;
  input Integer nc;
  input Real[nr,nc] Min;
  input Real[nq] q;
  output Real[nr, nr] Mout;
  protected
   Real[nr,nr] aux;
   Integer c=1;
 algorithm
   Mout := zeros(nr,nr);
   for i in 1:nq loop
     c :=(i - 1)*3;//column index
     aux[:,1] := Min[:,c+1];
     aux[:,2] := Min[:,c+2];
     aux[:,3] := Min[:,c+3];
     Mout := Mout+aux*q[i];
   end for;
 end getGrMatrix;

 function getGeMatrix
  input Integer nr;
  input Integer nq;
  input Integer nc;
  input Real[nr,nc] Min;
  input Real[nq] q;
  output Real[nq, 3] Mout;
  protected
   Real[nq,3] aux;
   Integer c=1;
 algorithm
   Mout := zeros(nq,3);
   for i in 1:nq loop
     c :=(i - 1)*3;//column index
     aux[:,1] := Min[:,c+1];
     aux[:,2] := Min[:,c+2];
     aux[:,3] := Min[:,c+3];
     Mout := Mout+aux*q[i];
   end for;
 end getGeMatrix;

  function getTaylorFunctionDerivative
   input Integer nr;
   input Integer nq;
   input Integer nc;
   input Real[nr,nc] M0;
   input Real[nr,nq,nc] M1;
   input Real[nq] q;
   output Real[nr,nc] M;
  protected
    Real[nr,nc] aux;
  algorithm
    M := zeros(nr,nc);
    for i in 1:nq loop
      aux := M1[:,i,:]*q[i];
      M := M+aux;
    end for;
    M := M0 +M;
   annotation(derivative=getTaylorFunctionDerivative);
  end getTaylorFunctionDerivative;
end MatrixFunctions;
