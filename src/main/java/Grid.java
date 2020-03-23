
import java.util.ArrayList;

public abstract class Grid
{
	ArrayList<Cell> cells;
	ArrayList<Face> faces;
	ArrayList<ScalarShapeFunction> shapeFunctions;
	DoubleTensor A;
	DoubleTensor rhs;
	public Grid()
	{
		cells = new ArrayList<>();
		faces = new ArrayList<>();
		shapeFunctions = new ArrayList<>();
	}
	public void evaluateCellIntegrals(ArrayList<CellIntegral> cellIntegrals,ArrayList<RightHandSideIntegral> rightHandSideIntegrals)
	{
		A = new DoubleTensor(shapeFunctions.size(),shapeFunctions.size(),true);
		rhs = new DoubleTensor(shapeFunctions.size());
		double integral;
		int i = 0;
		for(Cell K:cells)
		{
			System.out.println("evaluate cell integrals: "+(int)((1.0*i++)/(cells.size())*100)+"%");
			for(ScalarShapeFunction v:K.shapeFunctions)
			{
				for (ScalarShapeFunction u : K.shapeFunctions)
				{
					for(CellIntegral cellIntegral:cellIntegrals)
					{
						integral = cellIntegral.evaluateCellIntegral(K,u,v);
						A.add(u.globalIndex,v.globalIndex,integral);
					}

				}
				for(RightHandSideIntegral rightHandSideIntegral:rightHandSideIntegrals)
				{
					integral = rightHandSideIntegral.evaluateRightHandSideIntegral(K,v);
					rhs.add(v.globalIndex,integral);
				}
			}
		}
	}
	public void evaluateFaceIntegrals(ArrayList<FaceIntegral> faceIntegrals, ArrayList<BoundaryFaceIntegral> boundaryFaceIntegrals)
	{

		double integral;
		int i = 0;
		for(Face F:faces)
		{
			System.out.println("evaluate face integrals: "+(int)((1.0*i++)/(faces.size())*100)+"%");
			for(ScalarShapeFunction u:F.shapeFunctions)
			{
				for (ScalarShapeFunction v : F.shapeFunctions)
				{
					for(FaceIntegral faceIntegral:faceIntegrals)
					{
						integral = faceIntegral.evaluateFaceIntegral(F,u,v);
						A.add(u.globalIndex,v.globalIndex,integral);
					}

				}
			}
			if(F.isBoundaryFace)
			{
				for (ScalarShapeFunction v:F.shapeFunctions)
				{
					for(BoundaryFaceIntegral boundaryFaceIntegral:boundaryFaceIntegrals)
					{
						integral =
							boundaryFaceIntegral.evaluateBoundaryFaceIntegral(F, v);
						rhs.add(v.globalIndex, integral);
					}
				}
			}
		}
	}

}
