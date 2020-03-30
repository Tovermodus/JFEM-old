import com.google.common.collect.Multimap;

import java.util.ArrayList;
import java.util.stream.IntStream;

public class AngularGrid
{
	Grid g;
	DoubleTensor A;
	DoubleTensor TransAbs;
	DoubleTensor rhs;
	DoubleTensor Scat;
	DoubleTensor [] directions;
	double [] direction_weights;
	public AngularGrid(Grid grid, DoubleTensor[] directions, double [] direction_weights)
	{
		this.g = grid;
		this.directions = directions;
		this.direction_weights = direction_weights;
	}

	public AngularGrid refineAll(Multimap<Cell,Cell> cellRefinedCellMapping, Multimap<Face,Face> faceRefinedFaceMapping)
	{
		Grid refinedG = g.refineAll(cellRefinedCellMapping,faceRefinedFaceMapping);
		return new AngularGrid(refinedG,directions,direction_weights);
	}

	public void evaluateCellIntegrals(ArrayList<CellIntegral> cellIntegrals,
	                                  ArrayList<RightHandSideIntegral> rightHandSideIntegrals,
	                                  ArrayList<AngularCellIntegral> angularCellIntegrals)
	{
		A = new DoubleTensor(g.shapeFunctions.size()*directions.length,
			g.shapeFunctions.size()*directions.length,true);
		TransAbs = new DoubleTensor(g.shapeFunctions.size()*directions.length,
			g.shapeFunctions.size()*directions.length,true);
		Scat = new DoubleTensor(g.shapeFunctions.size()*directions.length,
			g.shapeFunctions.size()*directions.length,true);
		rhs = new DoubleTensor(g.shapeFunctions.size()*directions.length);
		g.evaluateCellIntegrals(cellIntegrals,rightHandSideIntegrals);
		IntStream.range(0,directions.length).parallel().forEach(i->
		{
			for(int j = 0; j < directions.length; j++)
			{
				for(Cell cell: g.cells)
				{
					for (ScalarShapeFunction u : cell.shapeFunctions)
					{
						for (ScalarShapeFunction v : cell.shapeFunctions)
						{
							double integral = 0;
							double nonScatIntegral = 0;
							for (AngularCellIntegral angularCellIntegral : angularCellIntegrals)
							{
								double inte = angularCellIntegral.evaluateAngularCellIntegral(cell,u,v,directions[i],directions[j],direction_weights[i],direction_weights[j]);
								integral += inte;
								if(!angularCellIntegral.name.equals(AngularCellIntegral.SCATTERING))
									nonScatIntegral += inte;
								else
									Scat.add(v.globalIndex*directions.length+j,
										u.globalIndex*directions.length+i,inte);
							}
							A.add(v.globalIndex*directions.length+j,
								u.globalIndex*directions.length+i,integral);
							if(nonScatIntegral != 0)
								TransAbs.add(v.globalIndex*directions.length+j,
									u.globalIndex*directions.length+i,nonScatIntegral);
						}
					}
				}
			}
		});

	}
	public void evaluateFaceIntegrals(ArrayList<FaceIntegral> faceIntegrals,
	                                  ArrayList<BoundaryFaceIntegral> boundaryFaceIntegrals,
	                                  ArrayList<AngularFaceIntegral> angularFaceIntegrals)
	{
		g.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		IntStream.range(0,directions.length).parallel().forEach(i->
		{
			for(int k = 0; k < g.A.sparseEntries; k++)
			{
				A.add(g.A.sparseYs[k]*directions.length+i,g.A.sparseXs[k]*directions.length+i,
					g.A.sparseValues[k]);
			}
			for(int k = 0; k < g.rhs.size(); k++)
			{
				rhs.add(k*directions.length+i,g.rhs.at(k));
			}
			for(int j = 0; j < directions.length; j++)
			{
				for(Face face: g.faces)
				{
					for (ScalarShapeFunction u : face.shapeFunctions)
					{
						for (ScalarShapeFunction v : face.shapeFunctions)
						{
							double integral = 0;
							for (AngularFaceIntegral angularFaceIntegral :
								angularFaceIntegrals)
							{
								integral += angularFaceIntegral.evaluateAngularFaceIntegral(face,u,v,
									directions[i],directions[j],direction_weights[i],direction_weights[j]);
							}
							A.add(v.globalIndex*directions.length+j,
								u.globalIndex*directions.length+i,integral);
							TransAbs.add(v.globalIndex*directions.length+j,
								u.globalIndex*directions.length+i,integral);
						}
					}
				}
			}
		});

	}
}
