within EMBSlib.Examples;

model Bracket
  Components.EMBS_Body eMBS_Body(SIDfileName = "modelica://EMBSlib/Resources/Data/bracket.SID_FEM") annotation(
    Placement(transformation(origin = {18, 0}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Mechanics.MultiBody.Parts.Fixed fixed annotation(
    Placement(transformation(origin = {-18, 0}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
  inner Modelica.Mechanics.MultiBody.World world annotation(
    Placement(transformation(origin = {-74, 40}, extent = {{-10, -10}, {10, 10}})));
equation
  connect(fixed.frame_b, eMBS_Body.frame_ref) annotation(
    Line(points = {{-8, 0}, {8, 0}}, color = {95, 95, 95}));
end Bracket;