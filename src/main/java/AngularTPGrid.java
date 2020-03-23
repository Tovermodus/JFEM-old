import java.util.ArrayList;

public class AngularTPGrid extends TPGrid
{
	DoubleTensor [] directions;
	double [] direction_weights;
	public AngularTPGrid(double xStart, double yStart, double xEnd, double yEnd, int numberXCells,
	                     int numberYCells, DoubleTensor [] directions, double [] direction_weights, int polynomialDegree)
	{
		super(xStart, yStart, xEnd, yEnd, numberXCells, numberYCells, polynomialDegree);
		this.directions = directions;
		this.direction_weights = direction_weights;
	}

	public void evaluateCellIntegrals(ArrayList<CellIntegral> cellIntegrals,
	                                  ArrayList<RightHandSideIntegral> rightHandSideIntegrals,
	                                  ArrayList<AngularTPCellIntegral> angularTPCellIntegrals)
	{
		A = new DoubleTensor(shapeFunctions.size()*directions.length,shapeFunctions.size()*directions.length,true);
		rhs = new DoubleTensor(shapeFunctions.size()*directions.length);
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
						for(int j = 0; j < directions.length; j++)
							A.add(u.globalIndex*directions.length + j,
								v.globalIndex*directions.length + j, integral);
					}

					for(AngularTPCellIntegral angularTPCellIntegral : angularTPCellIntegrals)
					{
						for(int j = 0; j < directions.length; j++)
							for(int k = 0; k < directions.length; k++)
							{
								integral =
									angularTPCellIntegral.evaluateAngularCellIntegral((TPCell) K, u, v, directions[j], directions[k], direction_weights[j], direction_weights[k]);
								A.add(u.globalIndex*directions.length+j,
									v.globalIndex*directions.length+k, integral);
							}

					}

				}
				for(RightHandSideIntegral rightHandSideIntegral:rightHandSideIntegrals)
				{
					integral =
						rightHandSideIntegral.evaluateRightHandSideIntegral(K,v);
					for(int j = 0; j < directions.length; j++)
						rhs.add(v.globalIndex*directions.length+j,integral);
				}

			}
		}
	}

	public void evaluateFaceIntegrals(ArrayList<FaceIntegral> faceIntegrals,
	                                  ArrayList<BoundaryFaceIntegral> boundaryFaceIntegrals,
	                                  ArrayList<AngularTPFaceIntegral> angularTPFaceIntegrals)
	{
		double value = 0;
		double integral = 0;
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
						for(int j = 0; j < directions.length; j++)
							A.add(u.globalIndex*directions.length+j,
								v.globalIndex*directions.length+j,
								integral);
					}
					for(AngularTPFaceIntegral angularFaceIntegral:angularTPFaceIntegrals)
					{
						for(int j = 0; j < directions.length; j++)
							for(int k = 0; k < directions.length; k++)
							{
								integral =
									angularFaceIntegral.evaluateAngularFaceIntegral((TPFace) F, u, v, directions[j],
								directions[k], direction_weights[j], direction_weights[k]);
								A.add(u.globalIndex
									*directions.length+j,
									v.globalIndex*directions.length+k,integral);
							}

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
						for(int j = 0; j < directions.length; j++)
							rhs.add(v.globalIndex*directions.length+j, integral);
					}
				}
			}
		}
	}
}
