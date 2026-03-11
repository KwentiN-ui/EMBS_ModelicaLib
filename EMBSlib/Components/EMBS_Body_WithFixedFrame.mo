within EMBSlib.Components;
model EMBS_Body_WithFixedFrame
  parameter Integer numNodes=9 "the number of nodes of the body";
  parameter Integer numModes=3 "the number of modes given in the SID file";
  parameter Real d=1 "damping scaling factor";
  parameter String SIDfileName=Modelica.Utilities.Files.loadResource("modelica://EMBSlib/Resources/Data/cartopPragV32.SID_FEM");
  parameter Boolean animation=true annotation(Dialog(tab="Animation"));
  parameter Real deformationScalingFactor=1 annotation(Dialog(tab="Animation"));
  parameter Real sphereDiameter=0.1 annotation(Dialog(tab="Animation"));
  parameter Integer[3] sphereColor={255,0,0} annotation(Dialog(tab="Animation"));
  parameter Real coordinateSystemScalingFactor=0.2 annotation(Dialog(tab="Animation"));
  constant Integer nr0=3 "dimension in space";
  Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_ref annotation(Placement(
   transformation(extent={{-116,-16},{-84,16}}),
   iconTransformation(extent={{-116,-16},{-84,16}})));
  Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_node[numNodes] annotation(Placement(
   transformation(extent={{84,-16},{116,16}}),
   iconTransformation(extent={{84,-16},{116,16}})));
  EMBSlib.Components.Internal.Node nodes[numNodes](
   each sid=sid,
   each nq=nq,
   nodeArrayIdx=1:numNodes,
   each deformationScalingFactor=deformationScalingFactor,
   each sphereColor=sphereColor,
   each sphereDiameter=sphereDiameter) annotation(Placement(transformation(extent={{-8,-4},{12,16}})));
  Real q[nq](start=zeros(nq)) "modal coordinates";
  Real qd[nq](start=zeros(nq))=der(q);
  Real qdd[nq](start=zeros(nq))=der(qd);
  Modelica.Units.SI.Position r_0[3](start={0,0,0})=frame_ref.r_0
        "position of frame of reference";
  Modelica.Units.SI.Velocity v[3]=der(r_0) "velocity of frame of reference";
  Modelica.Units.SI.Velocity v_0[3]=Modelica.Mechanics.MultiBody.Frames.resolve2(frame_ref.R,der(r_0))
        "velocity of frame of reference solved in body reference frame";
  Modelica.Units.SI.Acceleration a[3]=Modelica.Mechanics.MultiBody.Frames.resolve2(frame_ref.R,der(v))
        "acceleration of frame of reference solved in body reference frame";
  Modelica.Units.SI.Acceleration g_0_aux[3]=world.gravityAcceleration(frame_ref.r_0 + Modelica.Mechanics.MultiBody.Frames.resolve1(frame_ref.R,cm[:,1]));
  Modelica.Units.SI.Acceleration g_0[3]=Modelica.Mechanics.MultiBody.Frames.resolve2(frame_ref.R,g_0_aux)
        "Gravity acceleration resolved in body reference frame";
  Modelica.Units.SI.AngularVelocity omega[3]=Modelica.Mechanics.MultiBody.Frames.angularVelocity2(frame_ref.R);
  Modelica.Units.SI.AngularAcceleration omega_d[3]=der(omega);
  Modelica.Units.SI.AngularVelocity omega_tilde[3,3]={{0, -omega[3],omega[2]},  {omega[3],0,-omega[1]}, {-omega[2],omega[1],0}};
  Real M_t[3]=mI*(a-g_0) + transpose(mdCM_tilde)* omega_d + transpose(Ct)* qdd;
  Real k_omega_t[3] = mI*omega_tilde*v_0 + res2_1[:,1];
  Real res2_1[3,1]=omega_tilde*omega_tilde*mdCM; //the derivation of the center of gravity is neglected
  Real hd_t[3]=sum(identity(3)*nodes[i].f for i in 1:numNodes);
  Real M_r[3]=mdCM_tilde*(a-g_0) + J* omega_d + transpose(Cr) * qdd;
  Real k_omega_r[3]=Gr* omega + omega_tilde*J*omega + mdCM_tilde * omega_tilde * v_0;
  Real hd_r[3]=sum(identity(3)*nodes[i].t for i in 1:numNodes);
  Real M_q[nq]=Ct*(a-g_0) + Cr* omega_d + Me* qdd;
  Real k_omega_q[nq]=Ge*omega + Oe_*OMega + Ct * omega_tilde * v_0;
  Real k_q[nq]=ksigma[:,1] + Ke*q +d.*De*qd;
  Real hd_e[nq]=sum(nodes[i].hde_i for i in 1:numNodes);

  Modelica.Blocks.Sources.RealExpression qExp[nq](y=q) annotation(Placement(transformation(extent={{-58,40},{-38,60}})));
  Modelica.Units.SI.Torque[3] t_rest;
  Modelica.Units.SI.Force[3] f_rest;
  Modelica.Mechanics.MultiBody.Forces.WorldForce force(
   resolveInFrame=Modelica.Mechanics.MultiBody.Types.ResolveInFrameB.frame_b)
              annotation(Placement(transformation(extent={{-86,18},{-66,38}})));
  Modelica.Blocks.Sources.RealExpression f_elast1[nr0](y=f_rest) annotation(Placement(transformation(extent={{-124,16},{-104,36}})));
  Modelica.Mechanics.MultiBody.Forces.WorldTorque torque(
               resolveInFrame=Modelica.Mechanics.MultiBody.Types.ResolveInFrameB.world)
               annotation(Placement(transformation(extent={{-88,46},{-68,66}})));
  Modelica.Blocks.Sources.RealExpression t_elast2[nr0](y=t_rest) annotation(Placement(transformation(extent={{-124,46},{-104,66}})));
  protected
   parameter EMBSlib.Internal.SID_File sid=EMBSlib.Internal.SID_File(SIDfileName) annotation(Evaluate=true);
   parameter Modelica.Units.SI.Mass mass=EMBSlib.Internal.ExternalFunctions.getMass(sid);
   parameter Modelica.Units.SI.Mass mI[nr0,nr0]=identity(nr0)*mass;
   parameter Real mdCM_M0[nr0,1]=EMBSlib.Internal.ExternalFunctions.getM0( sid, "mdCM", nr0, 1) annotation(Evaluate=true);
   parameter Real mdCM_M1[nr0,nq,1]=EMBSlib.Internal.ExternalFunctions.getM1(sid,   "mdCM", nr0, nq, 1) annotation(Evaluate=true);
   Real mdCM[nr0,1]=EMBSlib.Internal.MatrixFunctions.getTaylorFunction(nr0,nq,1,mdCM_M0,mdCM_M1,q);
   Real cm[nr0,1]=mdCM/mass "centre of mass";
   Real mdCM_tilde[nr0,nr0]={{0, -mdCM[3,1],mdCM[2,1]},  {mdCM[3,1],0,-mdCM[1,1]}, {-mdCM[2,1],mdCM[1,1],0}};
   constant Integer nJ=6 "J is always a vektor of size 6";
   parameter Real J_M0[nJ,1]=EMBSlib.Internal.ExternalFunctions.getM0( sid, "J", nJ, 1) annotation(Evaluate=true);
   parameter Real J_M1[nJ,nq,1]=EMBSlib.Internal.ExternalFunctions.getM1( sid, "J", nJ, nq, 1) annotation(Evaluate=true);
   Real J_[nJ,1]=EMBSlib.Internal.MatrixFunctions.getTaylorFunction(nJ,nq,1,J_M0,J_M1,q);
   Real J[nr0,nr0]={{J_[1,1],J_[4,1],J_[5,1]},{J_[4,1],J_[2,1],J_[6,1]},{J_[5,1],J_[6,1],J_[3,1]}};
   parameter Real Ct_M0[nq,nr0]=EMBSlib.Internal.ExternalFunctions.getM0( sid, "Ct", nq, nr0) annotation(Evaluate=true);
   parameter Real Ct_M1[nq,nq,nr0]=EMBSlib.Internal.ExternalFunctions.getM1( sid, "Ct", nq, nq, nr0) annotation(Evaluate=true);
   Real Ct[nq,nr0]=EMBSlib.Internal.MatrixFunctions.getTaylorFunction(nq,nq,nr0,Ct_M0,Ct_M1,q);
   parameter Real Cr_M0[nq,nr0]=EMBSlib.Internal.ExternalFunctions.getM0( sid, "Cr", nq, nr0) annotation(Evaluate=true);
   parameter Real Cr_M1[nq,nq,nr0]=EMBSlib.Internal.ExternalFunctions.getM1( sid, "Cr", nq,  nq, nr0) annotation(Evaluate=true);
   Real Cr[nq,nr0]=EMBSlib.Internal.MatrixFunctions.getTaylorFunction(nq,nq,nr0,Cr_M0,Cr_M1,q);
   parameter Real Gr_[nr0,nr0*nq]=EMBSlib.Internal.ExternalFunctions.getM0( sid, "Gr", nr0, nr0*nq) annotation(Evaluate=true);
   Real Gr[nr0,nr0]=EMBSlib.Internal.MatrixFunctions.getGrMatrix(nr0,nq,nr0*nq,Gr_,qd);
   parameter Real Ge_[nq,nr0*nq]=EMBSlib.Internal.ExternalFunctions.getM0( sid, "Ge", nq, nr0*nq) annotation(Evaluate=true);
   Real Ge[nq,nr0]=EMBSlib.Internal.MatrixFunctions.getGeMatrix(nq,nq,nr0*nq,Ge_,qd);
   parameter Real Me[nq,nq]=identity(nq) annotation(Evaluate=true);
   constant Integer nOe=6 "its always 6";
   parameter Real Oe_M0[nq,nOe]=EMBSlib.Internal.ExternalFunctions.getM0( sid, "Oe", nq, nOe) annotation(Evaluate=true);
   parameter Real Oe_M1[nq,nq,nOe]=EMBSlib.Internal.ExternalFunctions.getM1( sid, "Oe", nq,  nq, nOe) annotation(Evaluate=true);
   Real Oe_[nq,nOe]=EMBSlib.Internal.MatrixFunctions.getTaylorFunction(nq,nq,6,Oe_M0,Oe_M1,q);
   Real OMega[nOe]={omega[1]^2, omega[2]^2, omega[3]^2, omega[1]*omega[2], omega[2]*omega[3], omega[1]*omega[3]};
   parameter Real ksigma[nq,1]=zeros(nq,1) annotation(Evaluate=true);
   parameter Real Ke[nq,nq]=EMBSlib.Internal.ExternalFunctions.getM0( sid, "Ke", nq, nq) annotation(Evaluate=true);
   parameter Real De[nq,nq]=EMBSlib.Internal.ExternalFunctions.getM0( sid, "De", nq, nq) annotation(Evaluate=true);
   parameter Integer nq=numModes;
   outer Modelica.Mechanics.MultiBody.World world annotation(Placement(transformation(extent={{-100,80},{-80,100}})));
  public
   Modelica.Mechanics.MultiBody.Visualizers.FixedFrame fixedFrame[numNodes](each
       length=coordinateSystemScalingFactor)                                                 annotation(Placement(transformation(extent={{60,60},{80,80}})));
  public
   Modelica.Mechanics.MultiBody.Visualizers.FixedFrame fixedFrame1[numNodes](
          each length=0.2)
        annotation (Placement(transformation(extent={{-28,-46},{-8,-26}})));
 equation
      //kinematic equations
      M_t + k_omega_t = -f_rest;
      M_r + k_omega_r = -t_rest;
      M_q + k_omega_q + k_q = hd_e;

      for i in 1:numNodes loop
       connect(qExp.y, nodes[i].q)    annotation (Line(points={{-37,50},{-18,
             50},{-18,5.8},{-8.2,5.8}},                                                                   color={0,0,127}));
       connect(nodes[i].frame_b, frame_node[i]) annotation (Line(
           points={{12,0.2},{25,0.2},{25,0},{100,0}},
           color={95,95,95},
           thickness=0.5));
       connect(nodes[i].frame_a, frame_ref) annotation (Line(
       points={{-8,0.4},{-55,0.4},{-55,0},{-100,0}},
       color={95,95,95},
       thickness=0.5));
      end for;

   connect(force.frame_b, frame_ref) annotation (Line(
       points={{-66,28},{-64,28},{-64,0},{-100,0}},
       color={95,95,95},
       thickness=0.5));
   connect(f_elast1.y,force.force) annotation(Line(
    points={{-103,26},{-98,26},{-93,26},{-93,28},{-88,28}},
    color={0,0,127}));
   connect(t_elast2.y,torque.torque) annotation(Line(
    points={{-103,56},{-98,56},{-95,56},{-90,56}},
    color={0,0,127}));
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
      connect(torque.frame_b, force.frame_b) annotation (Line(
          points={{-68,56},{-68,28},{-66,28}},
          color={95,95,95},
          thickness=0.5,
          smooth=Smooth.None));
  annotation (
   experiment(
    StopTime=1,
    StartTime=0,
    Interval=0.002,
    Algorithm="Dassl"));
 end EMBS_Body_WithFixedFrame;
