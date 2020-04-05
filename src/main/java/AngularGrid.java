import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import org.w3c.dom.UserDataHandler;

import java.awt.*;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;

public class AngularGrid
{
	Grid g;
	DoubleTensor A;
	DoubleTensor Areshuffled;
	ArrayList<DoubleTensor> cellMatrices;
	DoubleTensor rhs;
	DoubleTensor Scat;
	DoubleTensor shuffle;
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
		Areshuffled = new DoubleTensor(g.shapeFunctions.size()*directions.length,
			g.shapeFunctions.size()*directions.length,true);
		Scat = new DoubleTensor(g.shapeFunctions.size()*directions.length,
			g.shapeFunctions.size()*directions.length,true);
		rhs = new DoubleTensor(g.shapeFunctions.size()*directions.length);
		g.evaluateCellIntegrals(cellIntegrals,rightHandSideIntegrals);
		IntStream.range(0,directions.length).parallel().forEach(i->
		{

			List<List<Cell>> smallerList = Lists.partition(g.cells,12);
			smallerList.stream().parallel().forEach(smallList->
			{
				for (Cell cell : smallList)
				{
					for (int j = 0; j < directions.length; j++)
					{
						for (ScalarShapeFunction u : cell.shapeFunctions)
						{
							for (ScalarShapeFunction v : cell.shapeFunctions)
							{
								double integral = 0;
								double nonScatIntegral = 0;
								for (AngularCellIntegral angularCellIntegral : angularCellIntegrals)
								{
									double inte = angularCellIntegral.evaluateAngularCellIntegral(cell, u, v, directions[i], directions[j], direction_weights[i], direction_weights[j]);
									integral += inte;
									if (!angularCellIntegral.name.equals(AngularCellIntegral.SCATTERING))
										nonScatIntegral += inte;
									else
										Scat.add(v.globalIndex * directions.length + j,
											u.globalIndex * directions.length + i, inte);
								}
								if (integral != 0)
								{
									A.add(v.globalIndex * directions.length + j,
										u.globalIndex * directions.length + i, integral);
									Areshuffled.add(j*g.shapeFunctions.size()+v.globalIndex,
										i*g.shapeFunctions.size()+u.globalIndex, integral);
								}
							}
						}
					}
				}
			});
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

			List<List<Face>> smallerList = Lists.partition(g.faces,12);
			smallerList.stream().parallel().forEach(smallList->
			{
				for (Face face : smallList)
				{
					for (int j = 0; j < directions.length; j++)
					{
						for (ScalarShapeFunction u : face.shapeFunctions)
						{
							for (ScalarShapeFunction v : face.shapeFunctions)
							{
								double integral = 0;
								for (AngularFaceIntegral angularFaceIntegral :
									angularFaceIntegrals)
								{
									integral += angularFaceIntegral.evaluateAngularFaceIntegral(face, u, v,
										directions[i], directions[j], direction_weights[i], direction_weights[j]);
								}
								if (integral != 0)
								{
									A.add(v.globalIndex * directions.length + j,
										u.globalIndex * directions.length + i, integral);
									Areshuffled.add(j*g.shapeFunctions.size()+v.globalIndex,
										i*g.shapeFunctions.size()+u.globalIndex, integral);
								}
							}
						}
					}
				}
			});
		});
		shuffle = new DoubleTensor(A.getM(),A.getN(),true);
		for(int i = 0; i < g.shapeFunctions.size(); i++)
		{
			for(int j = 0; j < directions.length; j++)
			{
				shuffle.add(j*g.shapeFunctions.size()+i,i*directions.length+j,1);
			}
		}
	}
	public void generateCellMatrices()
	{
		g.cells.stream().parallel().forEach(cell ->
		{
			DoubleTensor cellMatrix = new DoubleTensor(cell.shapeFunctions.size()*directions.length,
				cell.shapeFunctions.size()*directions.length,false);
			for(int i = 0; i < cell.shapeFunctions.size(); i++)
			{
				for(int j = 0; j < cell.shapeFunctions.size(); j++)
				{
					for(int k = 0; k < directions.length; k++)
					{
						for(int l = 0; l < directions.length; l++)
						{
							ScalarShapeFunction u = cell.shapeFunctions.get(i);
							ScalarShapeFunction v = cell.shapeFunctions.get(j);

							cellMatrix.add(i*directions.length+k,j*directions.length+l,
								A.at(u.globalIndex*directions.length+k,
									v.globalIndex*directions.length+l));
						}
					}
				}
			}
			cellMatrices.add(cellMatrix);
		});
	}
}
