within EMBSlib.Components.Internal;
model Node
  parameter EMBSlib.Internal.SID_File sid;
  parameter Integer nq=3;
  parameter Integer nodeArrayIdx=1;
  parameter Integer nr0=3;
  parameter Real deformationScalingFactor=1;
  parameter Real coordinateSystemScalingFactor=1;
  parameter Integer[3] sphereColor={255,0,0};
  parameter Boolean animation=true annotation(Dialog(tab="Animation"));
  parameter Real origin_M0[nr0,1]=EMBSlib.Internal.ExternalFunctions.getM0Node( sid, "origin",nodeArrayIdx, nr0, 1) annotation(Evaluate=true);
  parameter Real origin_M1[nr0,nq,1]=EMBSlib.Internal.ExternalFunctions.getM1Node(sid, "origin",nodeArrayIdx, nr0, nq, 1) annotation(Evaluate=true);
  Real origin_[nr0,1]=origin_M0;
  Modelica.Units.SI.Position origin[nr0]=origin_[:,1];
  parameter Real psi_M0[nr0,nq]=EMBSlib.Internal.ExternalFunctions.getM0Node( sid, "psi",nodeArrayIdx, nr0, nq) annotation(Evaluate=true);
  parameter Real psi_M1[nr0,nq,nq]=EMBSlib.Internal.ExternalFunctions.getM1Node(sid, "psi",nodeArrayIdx, nr0, nq, nq) annotation(Evaluate=true);
  Real psi[nr0,nq]=psi_M0;
  Modelica.Units.SI.Angle theta[nr0]=psi*q "elastic rotation";
  Modelica.Units.SI.AngularVelocity der_theta[nr0]=der(theta)
        "elastic rotation velocity";
  parameter Real phi_M0[nr0,nq]=EMBSlib.Internal.ExternalFunctions.getM0Node( sid, "phi",nodeArrayIdx, nr0, nq) annotation(Evaluate=true);
  parameter Real phi_M1[nr0,nq,nq]=EMBSlib.Internal.ExternalFunctions.getM1Node(sid, "phi",nodeArrayIdx, nr0, nq, nq) annotation(Evaluate=true);
  Real phi[nr0,nq]=phi_M0;
  Modelica.Units.SI.Position u[nr0]=phi*q "elastic displacement";
  Modelica.Units.SI.Position u_abs=Modelica.Math.Vectors.length(u);
  parameter Real AP_M0[nr0,nr0]=EMBSlib.Internal.ExternalFunctions.getM0Node( sid, "AP", nodeArrayIdx, nr0, nr0) annotation(Evaluate=true);
  parameter Real AP_M1[nr0,nq,nr0]=EMBSlib.Internal.ExternalFunctions.getM1Node(sid, "AP",nodeArrayIdx, nr0, nq, nr0) annotation(Evaluate=true);
  Real AP[nr0,nr0]=EMBSlib.Internal.MatrixFunctions.getTaylorFunction(nr0,nq,nr0,AP_M0,zeros(nr0,nq,nr0),q);
  //Real AP[nr0,nr0]=EMBSlib.Internal.MatrixFunctions.getTaylorFunction(nr0,nq,nr0,AP_M0,AP_M1,q);
  Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a annotation(Placement(
   transformation(extent={{-116,-72},{-84,-40}}),
   iconTransformation(extent={{-116,-72},{-84,-40}})));
  Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b annotation(Placement(
   transformation(extent={{84,-74},{116,-42}}),
   iconTransformation(extent={{84,-74},{116,-42}})));
  Modelica.Blocks.Interfaces.RealInput q[nq](start=zeros(nq))
        "modal coordinates"                                                       annotation(Placement(
   transformation(extent={{-122,-22},{-82,18}}),
   iconTransformation(extent={{-122,-22},{-82,18}})));
  Real q_d[nq]=der(q);
  Modelica.Units.SI.Force f[nr0]=frame_b.f "external force applied";
  Modelica.Units.SI.Torque t[nr0]=frame_b.t "external torque applied";
  Modelica.Units.SI.Force hde_i[nq]=transpose(phi)*f+transpose(psi)*t;
  parameter Modelica.Units.SI.Diameter sphereDiameter=world.defaultBodyDiameter
        "Diameter of sphere"                                                                        annotation(Dialog(
   group="if animation = true",
   tab="Animation",
   enable=animation));
  Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape shape(
   shapeType="sphere",
   r=frame_b.r_0-sphereDiameter*0.5*{1,0,0},
   lengthDirection={1,0,0},
   length=sphereDiameter,
   width=sphereDiameter,
   height=sphereDiameter,
   color=sphereColor) annotation(Placement(transformation(extent={{80,60},{100,80}})));
  Modelica.Mechanics.MultiBody.Frames.Orientation R_theta=Modelica.Mechanics.MultiBody.Frames.axesRotations({1,2,3},theta,zeros(3));
 // Modelica.Mechanics.MultiBody.Frames.Orientation R_theta=Modelica.Mechanics.MultiBody.Frames.axesRotations({1,2,3},theta,der_theta); //unstable in some cases

  protected
   outer Modelica.Mechanics.MultiBody.World world;
 equation
   Connections.branch(frame_a.R, frame_b.R);
   assert(cardinality(frame_a) > 0 or cardinality(frame_b) > 0,
     "Neither connector frame_a nor frame_b of FixedTranslation object is connected");

   frame_b.r_0 = frame_a.r_0 + Modelica.Mechanics.MultiBody.Frames.resolve1(frame_a.R, origin+u);

   if Connections.rooted(frame_a.R) then
     frame_b.R = Modelica.Mechanics.MultiBody.Frames.absoluteRotation(frame_a.R, R_theta);
     zeros(3) = frame_a.f + Modelica.Mechanics.MultiBody.Frames.resolve1(R_theta, frame_b.f);
     //zeros(3) = frame_a.t + Modelica.Mechanics.MultiBody.Frames.resolve1(R_theta, frame_b.t) + cross(Modelica.Mechanics.MultiBody.Frames.resolve1(R_theta, origin+u),frame_b.f);
     zeros(3) = frame_a.t + Modelica.Mechanics.MultiBody.Frames.resolve1(R_theta, frame_b.t) - cross(origin+u,frame_a.f); //same result
   else
     frame_a.R = Modelica.Mechanics.MultiBody.Frames.absoluteRotation(frame_b.R, R_theta);
     zeros(3) = frame_b.f + Modelica.Mechanics.MultiBody.Frames.resolve1(R_theta, frame_a.f);
     zeros(3) = frame_b.t + Modelica.Mechanics.MultiBody.Frames.resolve1(R_theta, frame_a.t) + cross(Modelica.Mechanics.MultiBody.Frames.resolve1(R_theta, origin+u), frame_b.f);
   end if;
  annotation (
   Icon(
    coordinateSystem(preserveAspectRatio=false),
    graphics={
     Ellipse(
      fillColor={170,213,255},
      fillPattern=FillPattern.Solid,
      extent={{-54,50},{54,-58}}),
     Line(
      points={{36,36},{82,88}},
      thickness=0.5),
     Line(
      points={{-78,86},{-34,38}},
      thickness=0.5),
     Line(
      points={{-38,-44},{-86,-94}},
      thickness=0.5),
     Line(
      points={{32,-48},{94,-96}},
      thickness=0.5)}),
   Diagram(coordinateSystem(preserveAspectRatio=false)),
   experiment(
    StopTime=1,
    StartTime=0,
    Interval=0.001));
 end Node;
