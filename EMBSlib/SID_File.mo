within EMBSlib;
class SID_File
 extends ExternalObject;
 function constructor
  input String fileName;
  output SID_File sid;

  external "C" sid=SIDFileConstructor_C(fileName) 
    annotation(Library="modelica_rust",
               LibraryDirectory="modelica://EMBSlib/Resources/Library");
 end constructor;

 function destructor
  input SID_File sid;

  external "C" SIDFileDestructor_C(sid) 
    annotation(Library="modelica_rust",
               LibraryDirectory="modelica://EMBSlib/Resources/Library");
 end destructor;
end SID_File;
