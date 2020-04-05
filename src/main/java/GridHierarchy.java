import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Iterables;
import com.google.common.collect.Multimap;

import java.lang.reflect.Constructor;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.Stream;

public class GridHierarchy
{
	ArrayList<Grid> grids;
	ArrayList<Multimap<Cell,Cell>> cellMaps;
	ArrayList<Multimap<Face,Face>> faceMaps;
	ArrayList<DoubleTensor> prolongationMatrices;
	public GridHierarchy(Grid grid)
	{
		grids = new ArrayList<>();
		cellMaps = new ArrayList<>();
		faceMaps = new ArrayList<>();
		prolongationMatrices = new ArrayList<>();
		grids.add(grid);
	}
	public void addGloballyRefinedLevel()
	{
		Grid lastGrid = Iterables.getLast(grids);
		Multimap<Cell,Cell> cellMap = ArrayListMultimap.create();
		Multimap<Face,Face> faceMap = ArrayListMultimap.create();
		Grid refinedGrid = lastGrid.refineAll(cellMap,faceMap);
		grids.add(refinedGrid);
		cellMaps.add(cellMap);
		faceMaps.add(faceMap);
		DoubleTensor prolongationMatrix = new DoubleTensor(refinedGrid.shapeFunctions.size(),
			lastGrid.shapeFunctions.size(),true);
		prolongationMatrices.add(prolongationMatrix);
		for(ScalarShapeFunction shapeFunction:lastGrid.shapeFunctions)
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
				prolongationMatrix.add(i,shapeFunction.globalIndex,coeffs.get(i));
			}
		}
	}
	public void evaluateCellIntegrals(ArrayList<CellIntegral> cellIntegrals,ArrayList<RightHandSideIntegral> rightHandSideIntegrals)
	{
		for(Grid g:grids)
			g.evaluateCellIntegrals(cellIntegrals,rightHandSideIntegrals);
	}
	public void evaluateFaceIntegrals(ArrayList<FaceIntegral> faceIntegrals, ArrayList<BoundaryFaceIntegral> boundaryFaceIntegrals)
	{
		for(Grid g:grids)
			g.evaluateFaceIntegrals(faceIntegrals,boundaryFaceIntegrals);
	}

	public DoubleTensor multiGridSolve(double tolerance, int maxiter,
	                                                       Class<? extends Smoother>  smoother,String[] args)
	{
		ArrayList<Smoother> smoothers = new ArrayList<>();
		Constructor<? extends Smoother> constr = null;
		try
		{
			constr = smoother.getConstructor(Grid.class, String[].class);
		} catch (NoSuchMethodException e)
		{
			e.printStackTrace();
		}
		for (Grid grid : grids)
		{
			try
			{
				assert constr != null;
				smoothers.add(constr.newInstance(grid, args));
				System.out.println("kjh");
			} catch (Exception e)
			{
				e.printStackTrace();
			}
		}
		Grid lastGrid = Iterables.getLast(grids);
		System.out.println(lastGrid.A.getM()+" "+lastGrid.A.getN());
		DoubleTensor iterate = new DoubleTensor((int)lastGrid.rhs.size());
		DoubleTensor res = new DoubleTensor(lastGrid.rhs);
		DoubleTensor correct;
		System.out.println(res.vectorNorm());
		//res.print_formatted();
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
	public DoubleTensor multiGridV(int level, DoubleTensor rightHandSide, ArrayList<Smoother> smoothers)
	{
		int preIters = 2;
		int postIters = 2;
		DoubleTensor iterate = new DoubleTensor(rightHandSide.size());
		Grid g = grids.get(level);
		Smoother smoother = smoothers.get(level);
		//rightHandSide.print_formatted("rhs");
		if(level == 0)
		{
			//g.A.solve(rightHandSide).print_formatted("level0solve");
			return g.A.solve(rightHandSide);
		}
		for(int k = 0; k < preIters; k++)
			iterate = smoother.smooth(iterate,rightHandSide);
		DoubleTensor restrictedRightHandSide =
			prolongationMatrices.get(level - 1).tvmul(rightHandSide.sub(g.A.mvmul(iterate)));
		DoubleTensor correction = multiGridV(level - 1,
			restrictedRightHandSide,smoothers);
		iterate = iterate.add(prolongationMatrices.get(level - 1).mvmul(correction));
		//System.out.println("kajsdka");
		//System.out.println(rightHandSide.sub(g.A.mvmul(iterate)).vectorNorm());
		for(int k = 0; k < postIters; k++)
			iterate = smoother.smooth(iterate,rightHandSide);
		//System.out.println(rightHandSide.sub(g.A.mvmul(iterate)).vectorNorm());
		//iterate.print_formatted("return");
		return iterate;

	}

}
