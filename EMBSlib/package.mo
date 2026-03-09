within ;
package EMBSlib "An Elastic-Multi-Body Library"
  annotation(version="1.0.0",
             uses( Modelica(version="4.1.0")), Documentation(info="<html>
<p>The EMBSlib is an open-source Modelica library to model elastic bodies in a multi-body-system. The elastic bodies are parameterized by an SID-File which results from a reduced order finite element model.</p>
<p>Finite Element Models can be conducted in Ansys and exported as a reduced order model. With the help of the free PySID-tool or the commercial FEMBS-tool, SID-files can be generated which store the information on parameterize a flexible-body-model.</p>
<p>The elastic body can be used with the MultiBody-Library of the Modelica Standard Library.</p>
<p>Contributors: </p>
<p>Volker Waurich (TU-Dresden, Chair of Construction Machinery)</p>
<p>Robin Garbe (TU-Dresden, Chair of Construction Machinery)</p>
</html>"));
end EMBSlib;
