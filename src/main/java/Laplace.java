import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class Laplace
{
        public static void main(String[] args)
        {


                long startTime = System.nanoTime();

                System.out.println("output start");
                TPGrid grid = new TPGrid(0.,0.,1.,1.,1,1,1);
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
                System.out.println(((1.0*System.nanoTime() - startTime)/1e9));
                System.out.println("solve system: "+grid.A.getM()+"×"+grid.A.getN());
                //grid.A.makeParallelReady(12);
                DoubleTensor solution = grid.A.solveCG(grid.rhs,1e-11);
                System.out.println("solved");
                System.out.println(((1.0*System.nanoTime() - startTime)/1e9));
                //grid.A.print_formatted();
                //grid.rhs.print_formatted();
                FESpaceFunction solut = new FESpaceFunction(grid.shapeFunctions, solution);
                solut.plot(100,"/home/tovermodus/plot0.dat");
                System.out.println(((1.0*System.nanoTime() - startTime)/1e9));
        }
}
