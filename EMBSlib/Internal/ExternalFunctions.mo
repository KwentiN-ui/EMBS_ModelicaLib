within EMBSlib.Internal;
package ExternalFunctions
 function getMass
  input EMBSlib.Internal.SID_File sid;
  output Real mass;

  external "C" mass=getMass(sid) 
    annotation(Library="modelica_rust",
               LibraryDirectory="modelica://EMBSlib/Resources/Library");
 end getMass;

 function getNumNodes
  input EMBSlib.Internal.SID_File sid;
  output Integer numNodes;

  external "C" numNodes=getNumNodes(sid) 
    annotation(Library="modelica_rust",
               LibraryDirectory="modelica://EMBSlib/Resources/Library");
 end getNumNodes;

 function getNumModes
  input EMBSlib.Internal.SID_File sid;
  output Integer numModes;

  external "C" numModes=getNumModes(sid) 
    annotation(Library="modelica_rust",
               LibraryDirectory="modelica://EMBSlib/Resources/Library");
 end getNumModes;

 function getNumNodesFromFile
  input String fileName;
  output Integer numNodes;

  external "C" numNodes=getNumNodesFromFile(fileName) 
    annotation(Library="modelica_rust",
               LibraryDirectory="modelica://EMBSlib/Resources/Library");
 end getNumNodesFromFile;

 function getNumModesFromFile
  input String fileName;
  output Integer numModes;

  external "C" numModes=getNumModesFromFile(fileName) 
    annotation(Library="modelica_rust",
               LibraryDirectory="modelica://EMBSlib/Resources/Library");
 end getNumModesFromFile;

 function getM0
  input EMBSlib.Internal.SID_File sid;
  input String taylorName;
  input Integer nr;
  input Integer nc;
  output Real[nr,nc] m0;

  external "C" getM0(sid,taylorName,m0,nr,nc) 
    annotation(Library="modelica_rust",
               LibraryDirectory="modelica://EMBSlib/Resources/Library");
 end getM0;

 function getM1
  input EMBSlib.Internal.SID_File sid;
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
  input EMBSlib.Internal.SID_File sid;
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
  input EMBSlib.Internal.SID_File sid;
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
end ExternalFunctions;
