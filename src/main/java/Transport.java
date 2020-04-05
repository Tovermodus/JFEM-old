
import com.google.common.base.Stopwatch;
import com.google.common.collect.Iterables;

import java.security.DigestException;
import java.util.ArrayList;
import java.util.Arrays;

public class Transport
{
        public static void main(String[] args)
        {

//                DoubleTensor[] directions = new DoubleTensor[3];
//                double[] direction_weights = new double[directions.length];
//                for (int i = 0; i < directions.length; i++)
//                {
//                        double theta = 0.4+2. * i * Math.PI / directions.length;
//                        directions[i] = DoubleTensor.vectorFromValues(Math.sin(theta),
//                                Math.cos(theta));
//                        direction_weights[i] = 1;//2. * Math.PI / directions.length;
//                }
                DoubleTensor[] directions = Directions.TGLC2.getDirections();
                double[] direction_weights = Directions.TGLC2.getDirectionWeights();

                TPGrid grid = new TPGrid(-2, -2, 2, 2, 2, 2, 4);
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

                                return 0;
                        }
                };
                scatter.plot(100,directions,"/home/tovermodus/plots",-2,-2,2,2);
                AngularScalarFunction absorb = new AngularScalarFunction()
                {
                        @Override
                        public double value(DoubleTensor posInSpace, DoubleTensor posOnSphere1, DoubleTensor posOnSphere2)
                        {
                                return 0;
                        }
                };
                boolean diffLimit = false;
                angularCellIntegrals.add(new AngularTPCellIntegral(absorb,
                        AngularTPCellIntegral.ABSORPTION));
                angularCellIntegrals.add(new AngularTPCellIntegral(AngularTPCellIntegral.TRANSPORT));
                angularCellIntegrals.add(new AngularTPCellIntegral(scatter, AngularTPCellIntegral.SCATTERING));
                double penalty_param = 1.;


                if(diffLimit)
                {
                        angularFaceIntegrals.add(new AngularTPFaceIntegral(new AngularScalarFunction()
                        {
                                @Override
                                public double value(DoubleTensor posInSpace, DoubleTensor posOnSphere1, DoubleTensor
                                        posOnSphere2)
                                {

                                        return 4. / Math.max(4. / penalty_param, scatter.value(posInSpace, posOnSphere1,
                                                posOnSphere2) * (grid.xEnd - grid.xStart) / grid.numberXCells);
                                }
                        }, AngularTPFaceIntegral.JUMP_NORMALAVERAGE_JUMP_NORMALAVERAGE));
                        angularFaceIntegrals.add(new AngularTPFaceIntegral(AngularScalarFunction.constantFunction
                                (penalty_param),
                                AngularTPFaceIntegral.JUMP_NORMALAVERAGE_JUMP));
                        angularFaceIntegrals.add(new AngularTPFaceIntegral(AngularScalarFunction.constantFunction(penalty_param),
                                AngularTPFaceIntegral.VALUE_VALUE));
                }
                else
                {
                        angularFaceIntegrals.add(new AngularTPFaceIntegral(AngularScalarFunction.constantFunction(1)
                                ,AngularTPFaceIntegral.UPWIND));
                }
                rightHandSideIntegrals.add(new TPRightHandSideIntegral(new ScalarFunction()
                {
                        @Override
                        public double value(DoubleTensor pos)
                        {
                                double inner = 0.1;
                                double outer = 1;
                                double max = 200;
                                double min = 100;
                                if (pos.vectorNorm() < inner)
                                        return max+min;
                                else if (pos.vectorNorm() < outer)
                                        return min+(outer - pos.vectorNorm()) * max/(outer - inner);
                                return min;
                        }
                        @Override
                        public DoubleTensor derivative(DoubleTensor pos)
                        {
                                return null;
                        }
                }));
                gridHierarchy.addGloballyRefinedLevel();
                gridHierarchy.addGloballyRefinedLevel();
                gridHierarchy.addGloballyRefinedLevel();

                gridHierarchy.evaluateCellIntegrals(cellIntegrals,rightHandSideIntegrals,angularCellIntegrals);
                gridHierarchy.evaluateFaceIntegrals(faceIntegrals,boundaryFaceIntegrals,angularFaceIntegrals);
                AngularGrid ag = Iterables.getLast(gridHierarchy.angularGrids);
                System.out.println(ag.A.getM());
                System.out.println(ag.g.cells.size());
                int solMethod = 3;
                DoubleTensor solution2;
                LambdaSmoother gss = new LambdaSmoother(ag,new String[]{"0"});
                switch (solMethod)
                {
                        case 0:
                                ag.rhs.print_formatted();
                                solution2 = ag.A.solve(ag.rhs);
                                solution2.print_formatted();
                                break;
                        case 1:
                                solution2 = ag.A.solvePGMRES(new AngularMultigridPreconditioner(gridHierarchy,
                                gss.getClass(),new String[]{"0"}),ag.rhs,1e-5);
                                break;
                        case 2:

                                solution2 = gridHierarchy
                                        .multiGridSolve(1e-4,500,
                                                gss.getClass(),
                                                new String[]{"0"});
                                break;
                        case 3:

                                solution2 = new DoubleTensor(ag.rhs.size());
                                DoubleTensor res = new DoubleTensor(ag.rhs);
                                for(int i = 0; i < 30000&&res.vectorNorm()>1e-4;i++)
                                {
                                        res = ag.rhs.sub(ag.A.mvmul(solution2));
                                        System.out.println(i+" " +res.vectorNorm());
                                        solution2 = gss.smooth(solution2, ag.rhs);
                                }
                                break;
                        case 4:
                                solution2 = ag.A.solveGMRES(ag.rhs,1e-1);
                                break;
                        default:
                                throw new IllegalStateException("Unexpected value: " + solMethod);
                }
                AngularFESpaceFunction solfunc =
                        new AngularFESpaceFunction(Iterables.getLast(gridHierarchy.angularGrids).g.shapeFunctions,
                        directions,solution2);
                solfunc.plot(100,directions,"/home/tovermodus/plot",grid.xStart,grid.yStart,grid.xEnd,grid.yEnd);
                System.out.println("done");
        }
}
