within EMBSlib.Examples;

model RotatingBeam
  inner Modelica.Mechanics.MultiBody.World world annotation(
    Placement(transformation(origin = {-62, 0}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute revolute(n = {1, 0, 0}, useAxisFlange = true)  annotation(
    Placement(transformation(origin = {-10, 0}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Mechanics.Rotational.Sources.Speed speed(useSupport = true)  annotation(
    Placement(transformation(origin = {-16, 34}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Ramp ramp(duration = 1, height = 30000)  annotation(
    Placement(transformation(origin = {-58, 34}, extent = {{-10, -10}, {10, 10}})));
  Components.EMBS_Body eMBS_Body(SIDfileName = "modelica://EMBSlib/Resources/Data/beam_custom.SID_FEM")  annotation(
    Placement(transformation(origin = {46, 0}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Mechanics.MultiBody.Sensors.CutTorque cutTorque annotation(
    Placement(transformation(origin = {-36, 0}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Mechanics.MultiBody.Sensors.CutTorque cutTorque_lok annotation(
    Placement(transformation(origin = {18, 0}, extent = {{-10, -10}, {10, 10}})));
equation
  connect(speed.support, revolute.support) annotation(
    Line(points = {{-16, 24}, {-16, 10}}));
  connect(speed.flange, revolute.axis) annotation(
    Line(points = {{-6, 34}, {0, 34}, {0, 10}, {-10, 10}}));
  connect(ramp.y, speed.w_ref) annotation(
    Line(points = {{-46, 34}, {-28, 34}}, color = {0, 0, 127}));
  connect(world.frame_b, cutTorque.frame_a) annotation(
    Line(points = {{-52, 0}, {-46, 0}}, color = {95, 95, 95}));
  connect(cutTorque.frame_b, revolute.frame_a) annotation(
    Line(points = {{-26, 0}, {-20, 0}}, color = {95, 95, 95}));
  connect(revolute.frame_b, cutTorque_lok.frame_a) annotation(
    Line(points = {{0, 0}, {8, 0}}, color = {95, 95, 95}));
  connect(cutTorque_lok.frame_b, eMBS_Body.frame_ref) annotation(
    Line(points = {{28, 0}, {36, 0}}, color = {95, 95, 95}));
  annotation(
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.0001));
end RotatingBeam;