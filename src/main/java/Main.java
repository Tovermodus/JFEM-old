
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
                DoubleTensor[] directions = new DoubleTensor[4];
                double[] direction_weights = new double[directions.length];
                for (int i = 0; i < directions.length; i++)
                {
                        double theta = 2. * i * Math.PI / directions.length;
                        directions[i] = DoubleTensor.vectorFromValues(Math.sin(theta),
                                Math.cos(theta));
                        direction_weights[i] = 2. * Math.PI / directions.length;
                }
                AngularTPGrid grid = new AngularTPGrid(0, 0, 1, 1, 10, 10, directions, direction_weights, 2);
                ArrayList<CellIntegral> cellIntegrals = new ArrayList<>();
                ArrayList<FaceIntegral> faceIntegrals = new ArrayList<>();
                ArrayList<BoundaryFaceIntegral> boundaryFaceIntegrals = new ArrayList<>();
                ArrayList<RightHandSideIntegral> rightHandSideIntegrals = new ArrayList<>();
                ArrayList<AngularTPCellIntegral> angularTPCellIntegrals = new ArrayList<>();
                ArrayList<AngularTPFaceIntegral> angularTPFaceIntegrals = new ArrayList<>();
                AngularScalarFunction scatter = new AngularScalarFunction()
                {
                        @Override
                        public double value(DoubleTensor posInSpace, DoubleTensor posOnSphere1, DoubleTensor posOnSphere2)
                        {
                                return 1;
                        }
                };
                angularTPCellIntegrals.add(new AngularTPCellIntegral(AngularScalarFunction.constantFunction(5),
                        AngularTPCellIntegral.ABSORPTION));
                angularTPCellIntegrals.add(new AngularTPCellIntegral(AngularTPCellIntegral.TRANSPORT));
                angularTPCellIntegrals.add(new AngularTPCellIntegral(scatter, AngularTPCellIntegral.SCATTERING));
                double penalty_param = 10;
                angularTPFaceIntegrals.add(new AngularTPFaceIntegral(new AngularScalarFunction()
                {
                        @Override
                        public double value(DoubleTensor posInSpace, DoubleTensor posOnSphere1, DoubleTensor
                                posOnSphere2)
                        {

                                return penalty_param * 4. / Math.max(4., scatter.value(posInSpace, posOnSphere1,
                                        posOnSphere2)*(grid.xEnd - grid.xStart)/grid.numberXCells);
                        }
                }, AngularTPFaceIntegral.JUMP_NORMALAVERAGE_JUMP_NORMALAVERAGE));
                angularTPFaceIntegrals.add(new AngularTPFaceIntegral(AngularScalarFunction.constantFunction(penalty_param),
                        AngularTPFaceIntegral.JUMP_NORMALAVERAGE_JUMP));
                angularTPFaceIntegrals.add(new AngularTPFaceIntegral(AngularTPFaceIntegral.VALUE_VALUE));
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
                grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals, angularTPCellIntegrals);
                grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals, angularTPFaceIntegrals);
                System.out.println(grid.A.getM() + "×" + grid.A.getN());
                DoubleTensor solution = grid.A.transpose().solveGMRES(grid.rhs);
                System.out.println("solved");
                solution.print_formatted();
                int pointres = 100;
                double[][][] values = new double[pointres][pointres][directions.length];
                for (int k = 0; k < pointres; k++)
                {
                        for (int j = 0; j < pointres; j++)
                        {
                                for (int i = 0; i < grid.shapeFunctions.size(); i++)
                                {
                                        for (int l = 0; l < directions.length; l++)
                                        {
                                                if (i != -1)
                                                {
                                                        DoubleTensor pos =
                                                                DoubleTensor.vectorFromValues((grid.xEnd - grid.xStart) / (pointres - 1) * k,
                                                                        (grid.yEnd - grid.yStart) / (pointres - 1) * j);
                                                        values[j][k][l] += grid.shapeFunctions.get(i).value(pos) * solution.at(i * directions.length + l);
                                                }
                                        }
                                }
                        }
                }
                try
                {
                        for (int k = 0; k < directions.length; k++)
                        {
                                BufferedWriter plotWriter =
                                        new BufferedWriter(new FileWriter("/home/tovermodus/plot" + k + ".dat"));
                                for (int i = 0; i < pointres; i++)
                                {
                                        for (int j = 0; j < pointres; j++)
                                        {
                                                plotWriter.write(Double.toString(values[i][j][k]) + " ");
                                        }
                                        plotWriter.newLine();
                                }
                                plotWriter.flush();
                                plotWriter.close();
                        }
                } catch (IOException e)
                {
                        e.printStackTrace();
                }
                System.out.println("done");
        }
}
