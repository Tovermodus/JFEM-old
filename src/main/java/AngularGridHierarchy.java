import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Iterables;
import com.google.common.collect.Multimap;

import java.lang.reflect.Constructor;
import java.util.ArrayList;
import java.util.Map;

public class AngularGridHierarchy
{
	ArrayList<AngularGrid> angularGrids;
	ArrayList<Multimap<Cell,Cell>> cellMaps;
	ArrayList<Multimap<Face,Face>> faceMaps;
	ArrayList<DoubleTensor> prolongationMatrices;
	public AngularGridHierarchy(AngularGrid grid)
	{
		angularGrids = new ArrayList<>();
		cellMaps = new ArrayList<>();
		faceMaps = new ArrayList<>();
		prolongationMatrices = new ArrayList<>();
		angularGrids.add(grid);
	}
	public void addGloballyRefinedLevel()
	{
		AngularGrid lastGrid = Iterables.getLast(angularGrids);
		Multimap<Cell,Cell> cellMap = ArrayListMultimap.create();
		Multimap<Face,Face> faceMap = ArrayListMultimap.create();
		AngularGrid refinedGrid = lastGrid.refineAll(cellMap,faceMap);
		angularGrids.add(refinedGrid);
		cellMaps.add(cellMap);
		faceMaps.add(faceMap);
		DoubleTensor prolongationMatrix =
			new DoubleTensor(refinedGrid.g.shapeFunctions.size()*refinedGrid.directions.length,
			lastGrid.g.shapeFunctions.size()*lastGrid.directions.length,true);
		prolongationMatrices.add(prolongationMatrix);
		for(ScalarShapeFunction shapeFunction:lastGrid.g.shapeFunctions)
		{
			ArrayList<ScalarShapeFunction> refinedFunctions = new ArrayList<>();
			for(Cell cell:shapeFunction.cells)
			{
				for(Cell refinedSupportCell: cellMap.get(cell))
				{
					refinedFunctions.addAll(refinedSupportCell.shapeFunctions);
				}
			}
			Map<Integer, Double> coeffs = shapeFunction.prolongate(refinedFunctions);
			for(int i:coeffs.keySet())
			{
				for(int j = 0; j < lastGrid.directions.length; j++)
					prolongationMatrix.add(i*refinedGrid.directions.length+j,
						shapeFunction.globalIndex*lastGrid.directions.length+j,
						coeffs.get(i));
			}
		}
	}
	public void evaluateCellIntegrals(ArrayList<CellIntegral> cellIntegrals,
	                                  ArrayList<RightHandSideIntegral> rightHandSideIntegrals,
	                                  ArrayList<AngularCellIntegral> angularCellIntegrals)
	{
		for(AngularGrid g:angularGrids)
			g.evaluateCellIntegrals(cellIntegrals,rightHandSideIntegrals,angularCellIntegrals);
	}
	public void evaluateFaceIntegrals(ArrayList<FaceIntegral> faceIntegrals,
	                                  ArrayList<BoundaryFaceIntegral> boundaryFaceIntegrals,
	                                  ArrayList<AngularFaceIntegral> angularFaceIntegrals)
	{
		for(AngularGrid g:angularGrids)
			g.evaluateFaceIntegrals(faceIntegrals,boundaryFaceIntegrals,angularFaceIntegrals);
	}

	public DoubleTensor multiGridSolve(double tolerance, int maxiter,
	                                                       Class<? extends AngularSmoother>  smoother,String[] args)
	{
		ArrayList<AngularSmoother> smoothers = new ArrayList<>();
		Constructor<? extends AngularSmoother> constr = null;
		try
		{
			constr = smoother.getConstructor(AngularGrid.class, String[].class);
		} catch (NoSuchMethodException e)
		{
			e.printStackTrace();
		}
		for (AngularGrid grid : angularGrids)
		{
			try
			{
				assert constr != null;
				smoothers.add(constr.newInstance(grid, args));
			} catch (Exception e)
			{
				e.printStackTrace();
			}
		}
		AngularGrid lastGrid = Iterables.getLast(angularGrids);
		System.out.println(lastGrid.g.cells.size());
		DoubleTensor iterate = new DoubleTensor((int)lastGrid.rhs.size());
		DoubleTensor res = new DoubleTensor(lastGrid.rhs);
		DoubleTensor correct;
		System.out.println(res.vectorNorm());
		res.print_formatted();
		for (int i = 0; i < maxiter && res.vectorNorm() > tolerance; i++)
		{
			correct = multiGridV(prolongationMatrices.size(), res, smoothers);
			//correct.print_formatted("correct");
			iterate = iterate.add(correct);
			//iterate.print_formatted("mgiterate");
			res = lastGrid.rhs.sub(lastGrid.A.mvmul(iterate));
			//res.print_formatted("res");
			System.out.println(i + " " + res.vectorNorm());
		}
		return iterate;
	}
	public DoubleTensor multiGridV(int level, DoubleTensor rightHandSide, ArrayList<AngularSmoother> angularSmoothers)
	{
		int preIters = 5;
		int postIters = 5;
		DoubleTensor iterate = new DoubleTensor((int)rightHandSide.size());
		AngularGrid g = angularGrids.get(level);
		AngularSmoother angularSmoother = angularSmoothers.get(level);
		//rightHandSide.print_formatted("rhs");
		if(level == 0)
		{
			//g.A.solve(rightHandSide).print_formatted("level0solve");
			return g.A.solve(rightHandSide);
		}
		for(int k = 0; k < preIters; k++)
			iterate = angularSmoother.smooth(iterate,rightHandSide);
		DoubleTensor restrictedRightHandSide =
			prolongationMatrices.get(level - 1).tvmul(rightHandSide.sub(g.A.mvmul(iterate)));
		DoubleTensor correction = multiGridV(level - 1,
			restrictedRightHandSide,angularSmoothers);
		iterate = iterate.add(prolongationMatrices.get(level - 1).mvmul(correction));
		//System.out.println("kajsdka");
		//System.out.println(rightHandSide.sub(g.A.mvmul(iterate)).vectorNorm());
		for(int k = 0; k < postIters; k++)
			iterate = angularSmoother.smooth(iterate,rightHandSide);
		System.out.println(rightHandSide.sub(g.A.mvmul(iterate)).vectorNorm());
		//iterate.print_formatted("return");
		return iterate;

	}

}
