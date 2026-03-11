within EMBSlib.Examples;
model Beam "A Beam"
 extends Modelica.Icons.Example;
 EMBSlib.Components.EMBS_Body eMBS_Body(
     SIDfileName="modelica://EMBSlib/Resources/Data/Beam_3m.SID_FEM")                                                          annotation(Placement(transformation(extent={{-4,0},{
            16,20}})));
inner Modelica.Mechanics.MultiBody.World world(
 label2="z",
 g=9.81,
 n(each displayUnit="1")={0,0,-1},
     animateGravity=false)    annotation(Placement(transformation(extent={{-80,0},
             {-60,20}})));
Modelica.Mechanics.MultiBody.Joints.Revolute revolute1(
 useAxisFlange=true,
 n={0,1,0},
 w(fixed=false)) annotation(Placement(transformation(extent={{-38,0},{-18,20}})));
   Modelica.Mechanics.Rotational.Sources.Position position1(useSupport=true)
     annotation (Placement(transformation(extent={{-45,61},{-25,81}})));
   Modelica.Blocks.Sources.Ramp ramp(
     startTime=2,
     offset=-1.2,
     height=0.5,
     duration=2)
     annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
 EMBSlib.Components.EMBS_Body eMBS_Body2(
     SIDfileName="modelica://EMBSlib/Resources/Data/Beam_3m.SID_FEM")                                                          annotation(Placement(transformation(extent={{68,0},{
             88,20}})));
Modelica.Mechanics.MultiBody.Joints.Revolute revolute2(
 useAxisFlange=true,
 n={0,1,0},
 w(fixed=false)) annotation(Placement(transformation(extent={{36,0},{56,20}})));
   Modelica.Mechanics.Rotational.Sources.Position position2(useSupport=true)
     annotation (Placement(transformation(extent={{31,61},{51,81}})));
   Modelica.Blocks.Sources.Ramp ramp1(
     startTime=2,
     duration=2,
     height=0.3,
     offset=0) annotation (Placement(transformation(extent={{0,60},{20,80}})));
equation
connect(revolute1.frame_a,world.frame_b) annotation(Line(
 points={{-38,10},{-38,10},{-55,10},{-60,10}},
 color={95,95,95},
 thickness=0.0625));
connect(revolute1.frame_b,eMBS_Body.frame_ref) annotation(Line(
 points={{-18,10},{-4,10}},
 color={95,95,95},
 thickness=0.0625));
connect(position1.support,revolute1.support) annotation(Line(
 points={{-35,61},{-35,61},{-35,25},{-34,25},{-34,20}},
 thickness=0.0625));
connect(position1.flange,revolute1.axis) annotation(Line(
 points={{-25,71},{-20,71},{-20,45},{-30,45},{-30,20},{-28,20}},
 thickness=0.0625));
   connect(ramp.y, position1.phi_ref) annotation (Line(
       points={{-59,70},{-56,70},{-56,71},{-47,71}},
       color={0,0,127},
       smooth=Smooth.None));
connect(position2.support,revolute2.support) annotation(Line(
 points={{41,61},{43,61},{43,45},{40,45},{40,20}},
 thickness=0.0625));
connect(position2.flange,revolute2.axis) annotation(Line(
 points={{51,71},{56,71},{56,31},{46,31},{46,20}},
 thickness=0.0625));
   connect(ramp1.y, position2.phi_ref) annotation (Line(
       points={{21,70},{26,70},{26,71},{29,71}},
       color={0,0,127},
       smooth=Smooth.None));
   connect(revolute2.frame_b, eMBS_Body2.frame_ref) annotation (Line(
       points={{56,10},{68,10}},
       color={95,95,95},
       thickness=0.5,
       smooth=Smooth.None));
   connect(revolute2.frame_a, eMBS_Body.frame_node[1]) annotation (Line(
       points={{36,10},{28,10},{28,8.8},{16,8.8}},
       color={95,95,95},
       thickness=0.5,
       smooth=Smooth.None));
   annotation (
experiment(StopTime=5, __Dymola_Algorithm="Dassl"),
Diagram(
 coordinateSystem(
  preserveAspectRatio=false,
  extent={{-100,-100},{100,100}}),
graphics));
end Beam;
