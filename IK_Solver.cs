using System;
using System.Collections;
using System.Collections.Generic;

using Rhino;
using Rhino.Geometry;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;


using MathNet;
//using MathNet.Symbolics;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
  using Math = System.Math;
//using Exponential = MathNet.Symbolics.Exponential;
using Matrix = MathNet.Numerics.LinearAlgebra.Double.Matrix;




/// <summary>
/// This class will be instantiated on demand by the Script component.
/// </summary>
public class Script_Instance : GH_ScriptInstance
{
#region Utility functions
  /// <summary>Print a String to the [Out] Parameter of the Script component.</summary>
  /// <param name="text">String to print.</param>
  private void Print(string text) { /* Implementation hidden. */ }
  /// <summary>Print a formatted String to the [Out] Parameter of the Script component.</summary>
  /// <param name="format">String format.</param>
  /// <param name="args">Formatting parameters.</param>
  private void Print(string format, params object[] args) { /* Implementation hidden. */ }
  /// <summary>Print useful information about an object instance to the [Out] Parameter of the Script component. </summary>
  /// <param name="obj">Object instance to parse.</param>
  private void Reflect(object obj) { /* Implementation hidden. */ }
  /// <summary>Print the signatures of all the overloads of a specific method to the [Out] Parameter of the Script component. </summary>
  /// <param name="obj">Object instance to parse.</param>
  private void Reflect(object obj, string method_name) { /* Implementation hidden. */ }
#endregion

#region Members
  /// <summary>Gets the current Rhino document.</summary>
  private readonly RhinoDoc RhinoDocument;
  /// <summary>Gets the Grasshopper document that owns this script.</summary>
  private readonly GH_Document GrasshopperDocument;
  /// <summary>Gets the Grasshopper script component that owns this script.</summary>
  private readonly IGH_Component Component;
  /// <summary>
  /// Gets the current iteration count. The first call to RunScript() is associated with Iteration==0.
  /// Any subsequent call within the same solution will increment the Iteration count.
  /// </summary>
  private readonly int Iteration;
#endregion

  /// <summary>
  /// This procedure contains the user code. Input parameters are provided as regular arguments,
  /// Output parameters as ref arguments. You don't have to assign output parameters,
  /// they will have a default value.
  /// </summary>
  private void RunScript(Point3d Robot_Origin, DataTree<Mesh> Robot_Type, Plane TCP, List<double> DH_Param, ref object oMotorAngle, ref object oAxisPlanes, ref object oRobotMesh)
  {

    //================================Input Plane Reorientation==================================

    Point3d TCP_Origin = new Point3d(TCP.OriginX, TCP.OriginY, TCP.OriginZ);

    Vector3d YZ_Plane_Vec = new Vector3d(1.0, 0.0, 0.0);

    //define the initplane as TCP
    Plane initPlane = new Plane(TCP_Origin, YZ_Plane_Vec);
    Plane TCPPlane = new Plane(TCP_Origin, YZ_Plane_Vec);
    Transform rotTrans = Transform.Rotation(Math.PI * 0.5, TCPPlane.ZAxis, TCP_Origin);
    initPlane.Rotate(Math.PI * 0.5, TCPPlane.ZAxis, TCP_Origin);
    Transform testTrans = Transform.ChangeBasis(initPlane, TCP);

    double Roll = -Math.Atan2(testTrans.M10, testTrans.M00);
    double Pitch = Math.Atan2(-testTrans.M20, Math.Sqrt(Math.Pow(testTrans.M21, 2) + Math.Pow(testTrans.M22, 2)));
    double Yaw = -Math.Atan2(testTrans.M21, testTrans.M22);



    //================================Inverse Kinematics==================================

    double DH_d1 = DH_Param[0];
    double DH_a1 = DH_Param[1];
    double DH_a2 = DH_Param[2];
    double DH_a3 = DH_Param[3];
    double DH_d2 = DH_Param[4];
    double DH_dg = DH_Param[5];
    double DH_a4S = DH_Param[6];


    double q1 = new double();
    double q2 = new double();
    double q3 = new double();
    double q4 = new double();
    double q5 = new double();
    double q6 = new double();
    double d90 = Math.PI / 2;

    Matrix<double> T01 = pose(q1, 1, 1, DH_d1);
    Matrix<double> T12 = pose(q2 - d90, -d90, DH_a1, 0);
    Matrix<double> T23 = pose(q3, 0, DH_a2, 0);
    Matrix<double> T34 = pose(q4, -d90, DH_a4S, DH_d2);
    Matrix<double> T45 = pose(q5, d90, 0, 0);
    Matrix<double> T56 = pose(q6, -d90, 0, 0);
    Matrix<double> T6g = pose(0, 0, 0, DH_dg);

    Matrix<double> T03 = T01 * T12 * T23;
    Matrix<double> R03 = T03.RemoveColumn(3).RemoveRow(3);
    Matrix<double> R03T = R03.Transpose();

    Matrix<double> T36 = T34 * T45 * T56;
    Matrix<double> R36 = T36.RemoveColumn(3).RemoveRow(3);

    Matrix<double> Rgu = (rotz(Math.PI) * roty(-Math.PI / 2)).Transpose();
    Matrix<double> RguT = Rgu.Transpose();

    double px ;
    double py ;
    double pz ;

    px = TCP.OriginX;
    py = TCP.OriginY;
    pz = TCP.OriginZ;

    Point3d oPoint = new Point3d(px, py, pz);
    Plane oPlane = new Plane(oPoint, Vector3d.ZAxis);
    Transform t_roll = Transform.Rotation(Roll, oPlane.XAxis, oPlane.Origin);
    oPlane.Transform(t_roll);
    Transform t_pitch = Transform.Rotation(Pitch, oPlane.YAxis, oPlane.Origin);
    oPlane.Transform(t_pitch);
    Transform t_yaw = Transform.Rotation(Yaw, oPlane.ZAxis, oPlane.Origin);
    oPlane.Transform(t_yaw);

    Matrix<double> R0u_eval = rotz(Yaw) * roty(Pitch) * rotx(Roll);

    Matrix<double> R0g_eval = R0u_eval * RguT;

    Point3d gripperPoint = new Point3d(px, py, pz);

    Point3d wrist_center = get_wrist_center(gripperPoint, R0g_eval);


    List<double> q1q2q3_List = get_first_three_angles(wrist_center);

    T01 = pose(q1q2q3_List[0], 0, 0, DH_d1);
    T12 = pose(q1q2q3_List[1] - d90, -d90, DH_a1, 0);
    T23 = pose(q1q2q3_List[2], 0, DH_a2, 0);
    T03 = T01 * T12 * T23;
    R03 = T03.RemoveColumn(3).RemoveRow(3);
    R03T = R03.Transpose();

    Matrix<double> R36_eval = R03T * R0g_eval;

    List<double> q4q5q6_List = get_last_three_angles(R36_eval);
    List<double> q1q6_list = new List<double>();
    q1q6_list.AddRange(q1q2q3_List);
    q1q6_list.AddRange(q4q5q6_List);


    //================================Build Robot Body ==================================
    //====================robot body movement_Forward kinematics==================================

    //build origins for each Axis

    Point3d RobotOrigin = new Point3d(Robot_Origin.X, Robot_Origin.Y, Robot_Origin.Z);
    Vector3d DH_d1_Vec = new Vector3d(0, 0, DH_d1);
    Point3d Axis2_pre_origin = RobotOrigin + DH_d1_Vec;
    Vector3d DH_a1_Vec = new Vector3d(DH_a1, 0, 0);
    Point3d Axis2_origin = Axis2_pre_origin + DH_a1_Vec;
    Vector3d DH_a2_Vec = new Vector3d(0, 0, DH_a2);
    Point3d Axis3_origin = Axis2_origin + DH_a2_Vec;
    Vector3d DH_a3_Vec = new Vector3d(DH_a3, 0, 0);
    Point3d Axis4_origin = Axis3_origin + DH_a3_Vec;
    Vector3d DH_d2_Vec = new Vector3d(DH_d2, 0, 0);
    Point3d Axis5_origin = Axis3_origin + DH_d2_Vec;
    Vector3d DH_dg_Vec = new Vector3d(DH_dg, 0, 0);
    Point3d Axis6_origin = Axis5_origin + DH_dg_Vec;

    //build lines connnecting each axis origin
    Line A1_Ap2_Line = new Line(RobotOrigin, Axis2_pre_origin);
    Line Ap2_A2_Line = new Line(Axis2_pre_origin, Axis2_origin);
    Line A2_A3_Line = new Line(Axis2_origin, Axis3_origin);
    Line A3_A4_Line = new Line(Axis3_origin, Axis4_origin);
    Line A4_A5_Line = new Line(Axis4_origin, Axis5_origin);
    Line A5_A6_Line = new Line(Axis5_origin, Axis6_origin);

    PolyCurve robotLineBody = new PolyCurve();
    robotLineBody.Append(A1_Ap2_Line);
    robotLineBody.Append(Ap2_A2_Line);
    robotLineBody.Append(A2_A3_Line);
    robotLineBody.Append(A3_A4_Line);
    robotLineBody.Append(A4_A5_Line);
    robotLineBody.Append(A5_A6_Line);

    // rotate Axis 1
    Plane Robot_Origin_Plane = new Plane(RobotOrigin, Vector3d.ZAxis);//A1_origin plane
    Plane Robot_reOriented_Plane = new Plane(RobotOrigin, Vector3d.ZAxis);//A1_reoriented plane
    Transform A1_rot = Transform.Rotation(q1q2q3_List[0], Robot_Origin_Plane.Normal, RobotOrigin);
    robotLineBody.Transform(A1_rot);
    Robot_reOriented_Plane.Transform(A1_rot);
    Curve[] curveArray;
    curveArray = robotLineBody.Explode();

    // reconstruct polycurve for A2 rotation
    PolyCurve robotLineBody_rot_forA2 = new PolyCurve();
    for (int i = 2; i < curveArray.Length; i++)
    {
      robotLineBody_rot_forA2.Append(curveArray[i]);
    }

    // rotate Axis 2
    Plane Robot_A2_Origin_Plane = new Plane(curveArray[1].PointAtEnd, curveArray[2].PointAtEnd, curveArray[1].PointAtStart);//A2_origin plane
    Plane A2_reOriented_Plane = new Plane(curveArray[1].PointAtEnd, curveArray[2].PointAtEnd, curveArray[1].PointAtStart);//A2_reoriented plane
    Robot_A2_Origin_Plane.Flip();
    A2_reOriented_Plane.Flip();
    Transform rotate90DegreeOnA2 = Transform.Rotation(Math.PI * 0.5, Robot_A2_Origin_Plane.Normal, Robot_A2_Origin_Plane.Origin);
    Robot_A2_Origin_Plane.Transform(rotate90DegreeOnA2);
    A2_reOriented_Plane.Transform(rotate90DegreeOnA2);
    Transform A2_rot = Transform.Rotation(q1q2q3_List[1], Robot_A2_Origin_Plane.Normal, Robot_A2_Origin_Plane.Origin);
    robotLineBody_rot_forA2.Transform(A2_rot);
    A2_reOriented_Plane.Transform(A2_rot);
    // reconstruct polycurve for A3 rotation
    Curve[] curveArrayForA3;
    curveArrayForA3 = robotLineBody_rot_forA2.Explode();
    PolyCurve robotLineBody_rot_forA3 = new PolyCurve();
    for (int i = 1; i < curveArrayForA3.Length; i++)
    {
      robotLineBody_rot_forA3.Append(curveArrayForA3[i]);
    }
    // rotate Axis 3
    Plane Robot_A3_Origin_Plane = new Plane(curveArrayForA3[0].PointAtEnd, curveArrayForA3[0].PointAtStart, curveArrayForA3[1].PointAtEnd); // A3_origin plane
    Plane A3_reOriented_Plane = new Plane(curveArrayForA3[0].PointAtEnd, curveArrayForA3[0].PointAtStart, curveArrayForA3[1].PointAtEnd);//A3 reoriented plane
    Robot_A3_Origin_Plane.Flip();
    A3_reOriented_Plane.Flip();
    Transform rotate90DegreeOnA3 = Transform.Rotation(-Math.PI * 0.5, Robot_A3_Origin_Plane.Normal, Robot_A3_Origin_Plane.Origin);
    Robot_A3_Origin_Plane.Transform(rotate90DegreeOnA3);
    A3_reOriented_Plane.Transform(rotate90DegreeOnA3);
    Transform A3_rot = Transform.Rotation(q1q2q3_List[2], Robot_A3_Origin_Plane.Normal, Robot_A3_Origin_Plane.Origin);
    robotLineBody_rot_forA3.Transform(A3_rot);
    A3_reOriented_Plane.Transform(A3_rot);

    // there's a A4 shift need to repositioning
    Transform A4sTrans = Transform.Translation(Robot_A3_Origin_Plane.XAxis * DH_a4S);
    Point3d A4ShiftedOrigin = new Point3d();
    A4ShiftedOrigin = curveArrayForA3[1].PointAtEnd;
    A4ShiftedOrigin.Transform(A4sTrans);
    double A4sCstPointParam;
    curveArrayForA3[0].ClosestPoint(A4ShiftedOrigin, out A4sCstPointParam);
    Point3d TempA4point = curveArrayForA3[0].PointAt(A4sCstPointParam);
    Vector3d A4NormalVec = A4ShiftedOrigin - TempA4point;
    A4NormalVec.Unitize();

    //rotate Axis 4
    Plane A4_inital_Plane = new Plane(A4ShiftedOrigin, A4NormalVec);

    Transform planarProjected = Transform.PlanarProjection(A4_inital_Plane);
    Robot_A3_Origin_Plane.XAxis.Transform(planarProjected);
    Vector3d crossVecYAxis = Vector3d.CrossProduct(A4_inital_Plane.Normal, Robot_A3_Origin_Plane.XAxis);
    Plane Robot_A4_Origin_Plane = new Plane(A4_inital_Plane.Origin, Robot_A3_Origin_Plane.XAxis, crossVecYAxis);// A4_origin plane
    Plane A4_reOriented_Plane = new Plane(A4_inital_Plane.Origin, Robot_A3_Origin_Plane.XAxis, crossVecYAxis);//A4 reoriented plane
    Robot_A4_Origin_Plane.Transform(A3_rot);
    A4_reOriented_Plane.Transform(A3_rot);
    Transform A4_rot = Transform.Rotation(q4q5q6_List[0], A4_reOriented_Plane.Normal, A4_reOriented_Plane.Origin);
    A4_reOriented_Plane.Transform(A4_rot);

    //Construct Axis 5 origin plane
    Plane Robot_A5_Origin_Plane =new Plane(A4_reOriented_Plane.Origin, A4_reOriented_Plane.XAxis, A4_reOriented_Plane.YAxis);// A5_origin plane
    Plane A5_reOriented_Plane = new Plane(A4_reOriented_Plane.Origin, A4_reOriented_Plane.XAxis, A4_reOriented_Plane.YAxis);//A5 reoriented plane
    Transform rotate90DegreeforA5 = Transform.Rotation(Math.PI * 0.5, A4_reOriented_Plane.XAxis, A4_reOriented_Plane.Origin);
    Robot_A5_Origin_Plane.Transform(rotate90DegreeforA5);
    A5_reOriented_Plane.Transform(rotate90DegreeforA5);
    double A5distToA4 = DH_d2 - DH_a3;
    Transform A5PlanePosition = Transform.Translation(A4_reOriented_Plane.Normal * A5distToA4);
    Robot_A5_Origin_Plane.Transform(A5PlanePosition);
    A5_reOriented_Plane.Transform(A5PlanePosition);

    //Rotate Axis 5
    Transform A5_rot = Transform.Rotation(q4q5q6_List[1], Robot_A5_Origin_Plane.Normal, Robot_A5_Origin_Plane.Origin);
    A5_reOriented_Plane.Transform(A5_rot);

    // Construct Axis 6 origin plane
    Transform A6PlanePosition = Transform.Translation(A5_reOriented_Plane.YAxis * DH_dg);
    Plane Robot_A6_Origin_Plane = new Plane(A5_reOriented_Plane.Origin, A5_reOriented_Plane.XAxis, A5_reOriented_Plane.YAxis);// A6_origin plane
    Plane A6_reOriented_Plane = new Plane(A5_reOriented_Plane.Origin, A5_reOriented_Plane.XAxis, A5_reOriented_Plane.YAxis);//A6 reoriented plane
    Robot_A6_Origin_Plane.Transform(A6PlanePosition);
    A6_reOriented_Plane.Transform(A6PlanePosition);

    //Rotate Axis 6
    Transform rotate90DegreeforA6 = Transform.Rotation(Math.PI * 0.5, Robot_A6_Origin_Plane.XAxis, Robot_A6_Origin_Plane.Origin);
    Robot_A6_Origin_Plane.Transform(rotate90DegreeforA6);
    A6_reOriented_Plane.Transform(rotate90DegreeforA6);
    Robot_A6_Origin_Plane.Flip();
    A6_reOriented_Plane.Flip();
    Transform rotate90DegreeforA6_2 = Transform.Rotation(Math.PI * 0.5, A6_reOriented_Plane.Normal, A6_reOriented_Plane.Origin);
    Robot_A6_Origin_Plane.Transform(rotate90DegreeforA6_2);
    A6_reOriented_Plane.Transform(rotate90DegreeforA6_2);
    Transform A6_rot = Transform.Rotation(q4q5q6_List[2], Robot_A6_Origin_Plane.Normal, Robot_A6_Origin_Plane.Origin);
    A6_reOriented_Plane.Transform(A6_rot);

    ///////////////// mesh body building///////////////////
    // Orient mesh based on Axis 1 plane
    Transform OrientA1 = Transform.PlaneToPlane(Robot_Origin_Plane, Robot_reOriented_Plane);
    DataTree<Mesh > OrientA1_Meshes_Tree = new DataTree<Mesh>();

    OrientA1_Meshes_Tree = OrientMesh(Robot_Type, OrientA1);

    // Orient mesh based on Axis 2 plane
    Transform OrientA2 = Transform.PlaneToPlane(Robot_A2_Origin_Plane, A2_reOriented_Plane);
    DataTree<Mesh> OrientA2_Meshes_Tree = new DataTree<Mesh>();

    OrientA2_Meshes_Tree = OrientMesh(OrientA1_Meshes_Tree, OrientA2);

    // Orient mesh based on Axis 3 plane
    Transform OrientA3 = Transform.PlaneToPlane(Robot_A3_Origin_Plane, A3_reOriented_Plane);
    DataTree<Mesh> OrientA3_Meshes_Tree = new DataTree<Mesh>();

    OrientA3_Meshes_Tree = OrientMesh(OrientA2_Meshes_Tree, OrientA3);

    // Orient mesh based on Axis 4 plane
    Transform OrientA4 = Transform.PlaneToPlane(Robot_A4_Origin_Plane, A4_reOriented_Plane);
    DataTree<Mesh> OrientA4_Meshes_Tree = new DataTree<Mesh>();

    OrientA4_Meshes_Tree = OrientMesh(OrientA3_Meshes_Tree, OrientA4);

    // Orient mesh based on Axis 5 plane
    Transform OrientA5 = Transform.PlaneToPlane(Robot_A5_Origin_Plane, A5_reOriented_Plane);
    DataTree<Mesh> OrientA5_Meshes_Tree = new DataTree<Mesh>();

    OrientA5_Meshes_Tree = OrientMesh(OrientA4_Meshes_Tree, OrientA5);

    // Orient mesh based on Axis 6 plane
    Transform OrientA6 = Transform.PlaneToPlane(Robot_A6_Origin_Plane, A6_reOriented_Plane);
    DataTree<Mesh> OrientA6_Meshes_Tree = new DataTree<Mesh>();

    OrientA6_Meshes_Tree = OrientMesh(OrientA5_Meshes_Tree, OrientA6);


    //out put planes and Meshes

    List<Plane> RobotAxisReorientPlaneList = new List<Plane>();
    DataTree<Mesh> RobotMeshOutputTree = new DataTree<Mesh>();
    RobotMeshOutputTree.AddRange(Robot_Type.Branch(0));
    RobotMeshOutputTree.AddRange(OrientA1_Meshes_Tree.Branch(0));
    RobotMeshOutputTree.AddRange(OrientA2_Meshes_Tree.Branch(0));
    RobotMeshOutputTree.AddRange(OrientA3_Meshes_Tree.Branch(0));
    RobotMeshOutputTree.AddRange(OrientA4_Meshes_Tree.Branch(0));
    RobotMeshOutputTree.AddRange(OrientA5_Meshes_Tree.Branch(0));
    RobotMeshOutputTree.AddRange(OrientA6_Meshes_Tree.Branch(0));

    RobotAxisReorientPlaneList.Add(Robot_reOriented_Plane);
    RobotAxisReorientPlaneList.Add(A2_reOriented_Plane);
    RobotAxisReorientPlaneList.Add(A3_reOriented_Plane);
    RobotAxisReorientPlaneList.Add(A4_reOriented_Plane);
    RobotAxisReorientPlaneList.Add(A5_reOriented_Plane);
    RobotAxisReorientPlaneList.Add(A6_reOriented_Plane);





    oMotorAngle = q1q6_list;
    oAxisPlanes = RobotAxisReorientPlaneList; ;
    oRobotMesh = RobotMeshOutputTree;

  }

  // <Custom additional code> 



  public DataTree<Mesh> OrientMesh(DataTree<Mesh> inputMeshTree, Transform OrientMaxtrix)
  {
    DataTree<Mesh> returnTree = new DataTree<Mesh>();

    for (int i = 1; i < inputMeshTree.BranchCount; i++)
    {
      GH_Path pth = new GH_Path(i);
      if (inputMeshTree.Branch(i).Count > 1)
      {
        inputMeshTree.Branch(i)[0].Transform(OrientMaxtrix);
        inputMeshTree.Branch(i)[1].Transform(OrientMaxtrix);
        returnTree.AddRange(inputMeshTree.Branch(i), pth);
      }
      else
      {
        inputMeshTree.Branch(i)[0].Transform(OrientMaxtrix);
        returnTree.AddRange(inputMeshTree.Branch(i), pth);
      }
    }



    return returnTree;

  }

  public Point3d get_wrist_center(Point3d gripper_point, Matrix<double> R0g)
  {
    double xu = gripper_point.X;
    double yu = gripper_point.Y;
    double zu = gripper_point.Z;

    double nx = R0g.Row(0).At(2);
    double ny = R0g.Row(1).At(2);
    double nz = R0g.Row(2).At(2);

    double xw = xu - 0.23 * nx;
    double yw = yu - 0.23 * ny;
    double zw = zu - 0.23 * nz;

    Point3d wristCenter = new Point3d(xw, yw, zw);

    return wristCenter;

  }
  public Matrix<double> rotx(double q)
  {
    double sq = Math.Sin(q);
    double cq = Math.Cos(q);
    Matrix<double> r = DenseMatrix.OfArray(new double[]
      {
      {1,0, 0},
      {0, cq,-sq},
      {0, sq, cq}
      });
    return r;
  }

  public Matrix<double> roty(double q)
  {
    double sq = Math.Sin(q);
    double cq = Math.Cos(q);
    Matrix<double> r = DenseMatrix.OfArray(new double[]
      {
      {cq,0, sq},
      {0, 1,0},
      {-sq, 0,cq}
      });
    return r;
  }

  public Matrix<double> rotz(double q)
  {
    double sq = Math.Sin(q);
    double cq = Math.Cos(q);
    Matrix<double> r = DenseMatrix.OfArray(new double[]
      {
      {cq,-sq, 0},
      {sq, cq,0},
      {0, 0, 1},
      });
    return r;
  }

  public Matrix<double>pose (double theta, double alpha, double a, double d)
  {
    double r11 = Math.Cos(theta);
    double r12 = -Math.Sin(theta);
    double r23 = -Math.Sin(alpha);
    double r33 = Math.Cos(alpha);
    double r21 = Math.Sin(theta) * Math.Cos(alpha);
    double r22 = Math.Cos(theta) * Math.Cos(alpha);
    double r31 = Math.Sin(theta) * Math.Sin(alpha);
    double r32 = Math.Cos(theta) * Math.Sin(alpha);
    double y = -d * Math.Sin(alpha);
    double z = d * Math.Cos(alpha);

    Matrix<double> T = DenseMatrix.OfArray(new double[]
      {
      {r11,r12, 0.0, a},
      {r21, r22,r23, y},
      {r31, r32, r33, z},
      {0.0, 0.0, 0.0, 1.0}
      });

    return T;
  }

  public double get_hypotenuse(double a, double b)
  {
    return (Math.Sqrt(a * a + b * b));
  }

  public double get_cosine_law_angle(double a, double b, double c)
  {
    double cos_gamma = (a * a + b * b - c * c) / (2 * a * b);
    double sin_gamma = Math.Sqrt(1 - cos_gamma * cos_gamma);
    double gamma = Math.Atan2(sin_gamma, cos_gamma);
    return gamma;
  }

  public List<double> get_first_three_angles(Point3d wrist_center)
  {

    List<double > q1q2q3_List = new List<double>();

    double x = wrist_center.X;
    double y = wrist_center.Y;
    double z = wrist_center.Z;

    double a1 = 0.35;
    double a2 = 1.25;
    double a3 = -0.054;
    double d1 = 0.75;
    double d4 = 1.1;
    double l = get_hypotenuse(d4, -a3);
    double phi = Math.Atan2(d4, -a3);

    double x_prime = get_hypotenuse(x, y);
    double mx = x_prime - a1;
    double my = z - d1;
    double m = get_hypotenuse(mx, my);
    double alpha = Math.Atan2(my, mx);

    double gamma = get_cosine_law_angle(l, a2, m);
    double beta = get_cosine_law_angle(m, a2, l);

    double q1 = Math.Atan2(y, x);
    double q2 = Math.PI / 2 - beta - alpha;
    double q3 = -(gamma - phi);

    q1q2q3_List.Add(q1);
    q1q2q3_List.Add(q2);
    q1q2q3_List.Add(q3);

    return q1q2q3_List;
  }

  public List<double> get_last_three_angles(Matrix<double> R)
  {

    List<double> q4q5q6_List = new List<double>();

    double sin_q4 = R.Row(2).At(2);
    double cos_q4 = -R.Row(0).At(2);

    double sin_q5 = Math.Sqrt(Math.Pow(R.Row(0).At(2), 2) + Math.Pow(R.Row(2).At(2), 2));
    double cos_q5 = R.Row(1).At(2);

    double sin_q6 = -R.Row(1).At(1);
    double cos_q6 = R.Row(1).At(0);

    double q4 = Math.Atan2(sin_q4, cos_q4);
    double q5 = Math.Atan2(sin_q5, cos_q5);
    double q6 = Math.Atan2(sin_q6, cos_q6);

    q4q5q6_List.Add(q4);
    q4q5q6_List.Add(q5);
    q4q5q6_List.Add(q6);

    return q4q5q6_List;
  }
  public class IK_dictClass
  {
    public IK_dictClass() //inital the class
    {
    }
    // Define Modified DH Transformation matrix
    public Matrix<double> DH_Table(double alpha, double a, double d, double q)
    {
      Matrix <double> matrix = DenseMatrix.OfArray(new double[]
        {
        {Math.Cos(q), -Math.Sin(q), 0, a},
        {Math.Sin(q) * Math.Cos(alpha), Math.Cos(q) * Math.Cos(alpha), -Math.Sin(alpha), -Math.Sin(alpha) * d},
        {Math.Sin(q) * Math.Sin(alpha), Math.Cos(q) * Math.Sin(alpha), Math.Cos(alpha), Math.Cos(alpha) * d},
        {0, 0, 0, 1}
        });

      return matrix;
    }
    //Define RPY rotation matrices
    public Matrix<double> ROT_x(double r)
    {
      Matrix<double> matrix = DenseMatrix.OfArray(new double[]
        {
        {1,0,0},
        {0,Math.Cos(r),-Math.Sin(r) },
        {0,Math.Sin(r),Math.Cos(r) }
        });
      return matrix;
    }

    public Matrix<double> ROT_y(double p)
    {
      Matrix<double> matrix = DenseMatrix.OfArray(new double[]
        {
        {Math.Cos(p),0,Math.Sin(p)},
        {0,1,0},
        {-Math.Sin(p),0,Math.Cos(p) }
        });
      return matrix;
    }

    public Matrix<double> ROT_z(double y)
    {
      Matrix<double> matrix = DenseMatrix.OfArray(new double[]
        {
        {Math.Cos(y),-Math.Sin(y),0},
        {Math.Sin(y),Math.Cos(y),0},
        {0,0,1}
        });
      return matrix;
    }



  }

  // </Custom additional code> 
}
