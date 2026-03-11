within EMBSlib.Components;
model EMBS_Body
  extends EMBSlib.Components.Internal.PartialEMBSBody;
  annotation (
   Icon(coordinateSystem(preserveAspectRatio=false)),
   Diagram(coordinateSystem(preserveAspectRatio=false)),
   experiment(
    StopTime=1,
    StartTime=0,
    Interval=0.002,
    Algorithm="Dassl"));
end EMBS_Body;
