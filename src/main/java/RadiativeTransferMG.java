
import com.google.common.base.Stopwatch;
import com.google.common.collect.Iterables;

import java.util.ArrayList;

public class RadiativeTransferMG
{
        public static void main(String[] args)
        {
                DoubleTensor[] directions = new DoubleTensor[3];
                double[] direction_weights = new double[directions.length];
                for (int i = 0; i < directions.length; i++)
                {
                        double theta = 0.1+2. * i * Math.PI / directions.length;
                        directions[i] = DoubleTensor.vectorFromValues(Math.sin(theta),
                                Math.cos(theta));
                        direction_weights[i] = 1;//2. * Math.PI / directions.length;
                }
                TPGrid grid = new TPGrid(0, 0, 1, 1, 2, 2, 1);
                AngularGrid angularGrid = new AngularGrid(grid,directions,direction_weights);
                AngularGridHierarchy gridHierarchy = new AngularGridHierarchy(angularGrid);
                ArrayList<CellIntegral> cellIntegrals = new ArrayList<>();
                ArrayList<FaceIntegral> faceIntegrals = new ArrayList<>();
                ArrayList<BoundaryFaceIntegral> boundaryFaceIntegrals = new ArrayList<>();
                ArrayList<RightHandSideIntegral> rightHandSideIntegrals = new ArrayList<>();
                ArrayList<AngularCellIntegral> angularCellIntegrals = new ArrayList<>();
                ArrayList<AngularFaceIntegral> angularFaceIntegrals = new ArrayList<>();
                AngularScalarFunction scatter = new AngularScalarFunction()
                {
                        @Override
                        public double value(DoubleTensor posInSpace, DoubleTensor posOnSphere1, DoubleTensor posOnSphere2)
                        {
                                return 0.4;
                        }
                };
                angularCellIntegrals.add(new AngularTPCellIntegral(AngularScalarFunction.constantFunction(5),
                        AngularTPCellIntegral.ABSORPTION));
                angularCellIntegrals.add(new AngularTPCellIntegral(AngularTPCellIntegral.TRANSPORT));
                angularCellIntegrals.add(new AngularTPCellIntegral(scatter, AngularTPCellIntegral.SCATTERING));
                double penalty_param = 100.;
                angularFaceIntegrals.add(new AngularTPFaceIntegral(new AngularScalarFunction()
                {
                        @Override
                        public double value(DoubleTensor posInSpace, DoubleTensor posOnSphere1, DoubleTensor
                                posOnSphere2)
                        {

                                return 4./ Math.max(4., scatter.value(posInSpace, posOnSphere1,
                                        posOnSphere2)*(grid.xEnd - grid.xStart)/grid.numberXCells);
                        }
                }, AngularTPFaceIntegral.JUMP_NORMALAVERAGE_JUMP_NORMALAVERAGE));
                angularFaceIntegrals.add(new AngularTPFaceIntegral(AngularScalarFunction.constantFunction
                 (penalty_param),
                        AngularTPFaceIntegral.JUMP_NORMALAVERAGE_JUMP));
                angularFaceIntegrals.add(new AngularTPFaceIntegral(AngularTPFaceIntegral.VALUE_VALUE));
//                 angularFaceIntegrals.add(new AngularTPFaceIntegral(AngularScalarFunction.constantFunction(100),
//                        AngularTPFaceIntegral.UPWIND));
                rightHandSideIntegrals.add(new TPRightHandSideIntegral(new ScalarFunction()
                {
                        @Override
                        public double value(DoubleTensor pos)
                        {
                                if (pos.sub(DoubleTensor.vectorFromValues(0.5, 0.5)).vectorNorm() < 0.1)
                                        return 1;
                                else if (pos.sub(DoubleTensor.vectorFromValues(0.5, 0.5)).vectorNorm() < 0.35)
                                        return (0.35 - pos.sub(DoubleTensor.vectorFromValues(0.5, 0.5)).vectorNorm()) * 4;
                                else
                                        return 0;
                        }

                        @Override
                        public DoubleTensor derivative(DoubleTensor pos)
                        {
                                return null;
                        }
                }));
                gridHierarchy.addGloballyRefinedLevel();
                gridHierarchy.addGloballyRefinedLevel();
                gridHierarchy.evaluateCellIntegrals(cellIntegrals,rightHandSideIntegrals,angularCellIntegrals);
                gridHierarchy.evaluateFaceIntegrals(faceIntegrals,boundaryFaceIntegrals,angularFaceIntegrals);
                //gridHierarchy
                // .multiGridSolve(1e-3,5,
                //AngularGaussSeidelSmoother.class,
                //new String[0]);
                AngularGrid ag = gridHierarchy.angularGrids.get(2);
                LambdaSmoother ssm = new LambdaSmoother(ag,new String[]{"0.4"});
                DoubleTensor solution2 = new DoubleTensor(ag.A.getM());
                for(int i = 0; i < 25; i++)
                {
                        System.out.println("aksjdh");
                        solution2 = ssm.smooth(solution2, ag.rhs);
                        //solution2.print_formatted();
                        //solution.print_formatted();
                        DoubleTensor res = ag.rhs.sub(ag.A.mvmul(solution2));
                        System.out.println(res.vectorNorm());
                }
                AngularFESpaceFunction solfunc =
                        new AngularFESpaceFunction(Iterables.getLast(gridHierarchy.angularGrids).g.shapeFunctions,
                        directions,
                solution2);
                solfunc.plot(100,directions,"/home/tovermodus/plot");
                System.out.println("done");
        }
}
