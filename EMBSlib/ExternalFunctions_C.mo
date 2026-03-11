within EMBSlib;
package ExternalFunctions_C
 function getMass
  input EMBSlib.SID_File sid;
  output Real mass;

  external "C" mass=getMass(sid) 
    annotation(Library="modelica_rust",
               LibraryDirectory="modelica://EMBSlib/Resources/Library");
 end getMass;

 function getM0
  input EMBSlib.SID_File sid;
  input String taylorName;
  input Integer nr;
  input Integer nc;
  output Real[nr,nc] m0;

  external "C" getM0(sid,taylorName,m0,nr,nc) 
    annotation(Library="modelica_rust",
               LibraryDirectory="modelica://EMBSlib/Resources/Library");
 end getM0;

 function getM1
  input EMBSlib.SID_File sid;
  input String taylorName;
  input Integer nr;
  input Integer nq;
  input Integer nc;
  output Real[nr,nq,nc] m1;

  external "C" getM1(sid,taylorName,m1,nr,nq,nc) 
    annotation(Library="modelica_rust",
               LibraryDirectory="modelica://EMBSlib/Resources/Library");
 end getM1;

 function getM0Node
  input EMBSlib.SID_File sid;
  input String taylorName;
  input Integer nodeIdx;
  input Integer nr;
  input Integer nc;
  output Real[nr,nc] m0;

  external "C" getM0Node(sid,taylorName,nodeIdx,m0,nr,nc) 
    annotation(Library="modelica_rust",
               LibraryDirectory="modelica://EMBSlib/Resources/Library");
 end getM0Node;

 function getM1Node
  input EMBSlib.SID_File sid;
  input String taylorName;
  input Integer nodeIdx;
  input Integer nr;
  input Integer nq;
  input Integer nc;
  output Real[nr,nq,nc] m1;

  external "C" getM1Node(sid,taylorName,nodeIdx,m1,nr,nq,nc) 
    annotation(Library="modelica_rust",
               LibraryDirectory="modelica://EMBSlib/Resources/Library");
 end getM1Node;
end ExternalFunctions_C;
