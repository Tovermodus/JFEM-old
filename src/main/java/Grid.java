
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import org.checkerframework.checker.units.qual.A;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ForkJoinPool;

public class Grid
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

	public Grid refineAll(Multimap<Cell,Cell> cellRefinedCellMapping, Multimap<Face,Face> faceRefinedFaceMapping)
	{
		Grid ret = new Grid();
		ArrayList<Cell> refinedCells = ret.cells;
		ArrayList<Face> refinedFaces = ret.faces;
		ArrayList<ScalarShapeFunction> refinedShapeFunctions = ret.shapeFunctions;
		for(Cell cell:cells)
		{
			ArrayList<Cell> refCells = cell.refine(refinedFaces);
			for(Cell refCell: refCells)
			{
				refinedCells.add(refCell);
				cellRefinedCellMapping.put(cell, refCell);
			}
		}
		for(Face f:faces)
		{
			ArrayList<Face> refFaces = f.refine(cellRefinedCellMapping);
			for(Face refFace:refFaces)
			{
				refinedFaces.add(refFace);
				faceRefinedFaceMapping.put(f,refFace);
			}
		}
		/*for(Face f:refinedFaces)
		{
			f.cells = new ArrayList<>();
			for (Cell c : refinedCells)
			{
				c.faces = new ArrayList<>();
			}
		}
		for(Face f:refinedFaces)
		{
			for(Cell c:refinedCells)
			{
				if (c.isInCell(f.center()))
				{
					f.cells.add(c);
					c.faces.add(f);
				}
			}
		}*/
		for(Cell refinedCell: refinedCells)
		{
			refinedCell.distributeFunctions(refinedShapeFunctions);
		}
		return ret;

	}

	public void evaluateCellIntegrals(ArrayList<CellIntegral> cellIntegrals,ArrayList<RightHandSideIntegral> rightHandSideIntegrals)
	{
		A = new DoubleTensor(shapeFunctions.size(),shapeFunctions.size(),true);
		rhs = new DoubleTensor(shapeFunctions.size());

		List<List<Cell>> smallerList = Lists.partition(cells,12);
		int i = 0;
		ForkJoinPool pool = new ForkJoinPool(4);
		//pool.submit(()-> dd
		smallerList.stream().parallel().forEach(smallList->
		{
			for (Cell K : smallList)
			{
				//System.out.println("evaluate cell integrals: "+(int)((1.0*i++)/(cells.size())*100)+"%");
				for (ScalarShapeFunction v : K.shapeFunctions)
				{
					for (ScalarShapeFunction u : K.shapeFunctions)
					{
						double integral = 0;
						for (CellIntegral cellIntegral : cellIntegrals)
						{
							integral += cellIntegral.evaluateCellIntegral(K, u, v);
						}
						if(integral != 0)
							A.add(v.globalIndex, u.globalIndex, integral);
					}
					double integral = 0;
					for (RightHandSideIntegral rightHandSideIntegral : rightHandSideIntegrals)
					{
						integral += rightHandSideIntegral.evaluateRightHandSideIntegral(K, v);
					}
					if(integral != 0)
						rhs.add(v.globalIndex, integral);

				}
			}

		});
	}
	public void evaluateFaceIntegrals(ArrayList<FaceIntegral> faceIntegrals, ArrayList<BoundaryFaceIntegral> boundaryFaceIntegrals)
	{

		List<List<Face>> smallerList = Lists.partition(faces,12);
		int i = 0;
		ForkJoinPool pool = new ForkJoinPool(4);
		//pool.submit(()-> dd
		smallerList.stream().parallel().forEach(smallList->
		{
			for (Face F : smallList)
			{
				//System.out.println("evaluate face integrals: " + (int) ((1.0 * i++) / (faces.size()
				// ) * 100) + "%");
				for (ScalarShapeFunction u : F.shapeFunctions)
				{
					for (ScalarShapeFunction v : F.shapeFunctions)
					{

						double integral = 0;
						for (FaceIntegral faceIntegral : faceIntegrals)
						{
							integral += faceIntegral.evaluateFaceIntegral(F,u,v);
						}
						if(integral != 0)
							A.add(v.globalIndex, u.globalIndex, integral);
					}
				}
				if (F.isBoundaryFace)
				{
					for (ScalarShapeFunction v : F.shapeFunctions)
					{
						double integral = 0;
						for (BoundaryFaceIntegral boundaryFaceIntegral : boundaryFaceIntegrals)
						{
							integral +=
								boundaryFaceIntegral.evaluateBoundaryFaceIntegral(F, v);
						}
						if(integral != 0)
							rhs.add(v.globalIndex, integral);
					}
				}
			}
		});
	}

}
