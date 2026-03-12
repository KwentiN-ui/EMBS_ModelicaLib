within EMBSlib.Examples;

model BeamCCX
  inner Modelica.Mechanics.MultiBody.World world annotation(
    Placement(transformation(origin = {-54, 0}, extent = {{-10, -10}, {10, 10}})));
  Components.EMBS_Body_WithFixedFrame eMBS_Body_WithFixedFrame(SIDfileName = "modelica://EMBSlib/Resources/Data/beam_custom.SID_FEM")  annotation(
    Placement(transformation(origin = {-24, 0}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Mechanics.MultiBody.Forces.WorldForce force annotation(
    Placement(transformation(origin = {8, 0}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  Modelica.Blocks.Sources.Step step[3](height = {0, -100, 0}, startTime = {0, 0.1, 0})  annotation(
    Placement(transformation(origin = {46, 0}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
equation
  connect(force.frame_b, eMBS_Body_WithFixedFrame.frame_node[6]) annotation(
    Line(points = {{-2, 0}, {-14, 0}}, color = {95, 95, 95}));
  connect(step.y, force.force) annotation(
    Line(points = {{36, 0}, {20, 0}}, color = {0, 0, 127}));
  connect(world.frame_b, eMBS_Body_WithFixedFrame.frame_ref) annotation(
    Line(points = {{-44, 0}, {-34, 0}}, color = {95, 95, 95}));
annotation(
    experiment(StartTime = 0, StopTime = 0.2, Tolerance = 1e-06, Interval = 2e-05));
end BeamCCX;