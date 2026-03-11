within EMBSlib.Components;
model EMBS_Body_WithFixedFrame
  extends EMBSlib.Components.Internal.PartialEMBSBody;
  public
   Modelica.Mechanics.MultiBody.Visualizers.FixedFrame fixedFrame[numNodes](each
       length=coordinateSystemScalingFactor)                                                 annotation(Placement(transformation(extent={{60,60},{80,80}})));
  public
   Modelica.Mechanics.MultiBody.Visualizers.FixedFrame fixedFrame1[numNodes](
          each length=0.2)
        annotation (Placement(transformation(extent={{-28,-46},{-8,-26}})));
 equation
   connect(fixedFrame.frame_a, nodes.frame_b) annotation (Line(
       points={{60,70},{38,70},{38,0.2},{12,0.2}},
       color={95,95,95},
       thickness=0.5,
       smooth=Smooth.None));
      connect(fixedFrame1.frame_a, nodes.frame_a) annotation (Line(
          points={{-28,-36},{-38,-36},{-38,-14},{-8,-14},{-8,0.4}},
          color={95,95,95},
          thickness=0.5,
          smooth=Smooth.None));
  annotation (
   Icon(coordinateSystem(preserveAspectRatio=false)),
   Diagram(coordinateSystem(preserveAspectRatio=false)),
   experiment(
    StopTime=1,
    StartTime=0,
    Interval=0.002,
    Algorithm="Dassl"));
end EMBS_Body_WithFixedFrame;
