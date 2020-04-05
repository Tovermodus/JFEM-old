
import com.google.common.base.Stopwatch;
import jeigen.SparseMatrixLil;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class Main
{
        public static void main(String[] args)
        {
                DoubleTensor[] directions = new DoubleTensor[3];
                double[] direction_weights = new double[directions.length];
                for (int i = 0; i < directions.length; i++)
                {
                        double theta = 2. * i * Math.PI / directions.length;
                        directions[i] = DoubleTensor.vectorFromValues(Math.sin(theta),
                                Math.cos(theta));
                        direction_weights[i] = 1;//2. * Math.PI / directions.length;
                }
                TPGrid grid = new TPGrid(0, 0, 1, 1, 12, 12, 1);
                AngularGrid angularGrid = new AngularGrid(grid,directions,direction_weights);
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
                                return 01.4;
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
                angularFaceIntegrals.add(new AngularTPFaceIntegral(AngularScalarFunction.constantFunction(penalty_param),
                        AngularTPFaceIntegral.JUMP_NORMALAVERAGE_JUMP));
                angularFaceIntegrals.add(new AngularTPFaceIntegral(AngularTPFaceIntegral.VALUE_VALUE));
                //angularTPFaceIntegrals.add(new AngularTPFaceIntegral(AngularTPFaceIntegral.JUMP_NORMALAVERAGE_JUMP_NORMALAVERAGE));
                // angularTPFaceIntegrals.add(new AngularTPFaceIntegral(AngularScalarFunction.constantFunction(1000),
                //        AngularTPFaceIntegral.UPWIND));
                //angularTPFaceIntegrals.add(new AngularTPFaceIntegral(AngularScalarFunction.constantFunction(10),
                //        AngularTPFaceIntegral.UPWIND));
                rightHandSideIntegrals.add(new TPRightHandSideIntegral(new ScalarFunction()
                {
                        @Override
                        public double value(DoubleTensor pos)
                        {
                                if (pos.sub(DoubleTensor.vectorFromValues(0.5, 0.5)).vectorNorm() < 0.1)
                                        return 100;
                                else if (pos.sub(DoubleTensor.vectorFromValues(0.5, 0.5)).vectorNorm() < 0.35)
                                        return (0.35 - pos.sub(DoubleTensor.vectorFromValues(0.5, 0.5)).vectorNorm()) * 400;
                                else
                                        return 0;
                        }

                        @Override
                        public DoubleTensor derivative(DoubleTensor pos)
                        {
                                return null;
                        }
                }));
                angularGrid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals, angularCellIntegrals);
                System.out.println("cells done");
                angularGrid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals, angularFaceIntegrals);
                System.out.println("faces done");
                System.out.println(angularGrid.A.getM() + "×" + angularGrid.A.getN());
                Stopwatch stop = Stopwatch.createStarted();
                DoubleTensor solution = angularGrid.A.solve(angularGrid.rhs);
                System.out.println(stop.elapsed().toMillis());
                System.out.println("solved");
                solution.print_formatted();
                AngularFESpaceFunction solfunc = new AngularFESpaceFunction(grid.shapeFunctions,directions,solution);
                //solfunc.plot(100,directions,"/home/tovermodus/plot");
                System.out.println("done");
        }
}
