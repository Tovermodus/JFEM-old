
import jeigen.SparseMatrixLil;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class Laplace
{
        public static void main(String[] args)
        {
                long startTime = System.nanoTime();

                System.out.println("output start");
                TPGrid grid = new TPGrid(0.,0.,1.,1.,32,32,1);
                CellIntegral gg = new TPCellIntegral(TPCellIntegral.GRAD_GRAD);
                TPFaceIntegral jj = new TPFaceIntegral(ScalarFunction.constantFunction(1000.0),
                        TPFaceIntegral.VALUE_JUMP_VALUE_JUMP);
                ArrayList<CellIntegral> cellIntegrals = new ArrayList<>();
                cellIntegrals.add(gg);
                ArrayList<FaceIntegral> faceIntegrals = new ArrayList<>();
                faceIntegrals.add(jj);
                TPRightHandSideIntegral rightHandSideIntegral =
                        new TPRightHandSideIntegral(ScalarFunction.oneFunction());
                ArrayList<RightHandSideIntegral> rightHandSideIntegrals = new ArrayList<>();
                rightHandSideIntegrals.add(rightHandSideIntegral);
                ArrayList<BoundaryFaceIntegral> boundaryFaceIntegrals = new ArrayList<>();
                grid.evaluateCellIntegrals(cellIntegrals,rightHandSideIntegrals);
                grid.evaluateFaceIntegrals(faceIntegrals,boundaryFaceIntegrals);
                System.out.println("solve system: "+grid.A.getM()+"×"+grid.A.getN());
                DoubleTensor solution = grid.A.solveCG(grid.rhs,1e-7);
                System.out.println("solved");
                System.out.println(((1.0*System.nanoTime() - startTime)/1e9));
                int pointres = 100;
                double[][] values = new double[pointres][pointres];
                for(int k = 0; k < pointres; k++)
                {
                        for(int j = 0; j < pointres; j++)
                        {
                                for(int i = 0; i <grid.shapeFunctions.size();i++)
                                {
                                        if(i !=-4)
                                        {
                                                DoubleTensor pos = DoubleTensor.vectorFromValues(1.0 / (pointres-1) * k,
                                                        1.0 / (pointres-1) * j);
                                                values[k][j] += grid.shapeFunctions.get(i).value(pos)*solution.at(i);
                                        }
                                }
                        }
                }
                try
                {
                        BufferedWriter plotWriter = new BufferedWriter(new FileWriter("/home/tovermodus/plot.dat"));
                        for(int i = 0; i < pointres; i++)
                        {
                                for(int j = 0; j < pointres; j++)
                                {
                                        plotWriter.write(Double.toString(values[i][j])+" ");
                                }
                                plotWriter.newLine();
                        }
                        plotWriter.flush();
                        plotWriter.close();
                } catch (IOException e)
                {
                        e.printStackTrace();
                }
        }
}
