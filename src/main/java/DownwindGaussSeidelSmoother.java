import com.google.common.collect.BiMap;

import java.util.ArrayList;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class DownwindGaussSeidelSmoother extends  AngularSmoother
{

	BlockMatrix blockMatrix;
	BlockDiagonalMatrix blockDiag;
	BlockDiagonalMatrix blockInverse;
	ArrayList<Integer> downWindOrderings[];
	int blockSize;
	int blockN;
	public DoubleTensor getCellMatrix(int cellIndex1, int cellIndex2, BlockMatrix m)
	{
		return m.blocks[cellIndex1][cellIndex2];
	}
	public DoubleTensor getCellVector(int cellIndex, DoubleTensor vector)
	{
		DoubleTensor cellVector = new DoubleTensor(blockSize);
		for(int i = 0; i < blockSize; i++)
		{
			cellVector.set(i,vector.at(cellIndex*blockSize+i));
		}
		return cellVector;
	}
	public void setCellVector(int cellindex, DoubleTensor largeVector, DoubleTensor cellVector)
	{
		for(int i = 0; i < blockSize; i++)
		{
			largeVector.set(cellindex*blockSize+i,cellVector.at(i));
		}
	}
	public ArrayList<Integer> getDownwindOrdering(DoubleTensor dir)
	{
		ArrayList<Integer> ret = new ArrayList<>();
		for(int i = 0; i < grid.g.cells.size();i++)
			ret.add(i);
		ret.sort((i,j)->{double inner =
			grid.g.cells.get(i).center().sub(grid.g.cells.get(j).center()).inner(dir);
		if(inner > 1e-14)
		return 1;
		else if(inner < -1e-14)
		return -1;
		else
		return 0;});
		return ret;
	}
	public DownwindGaussSeidelSmoother(AngularGrid g, String[] args)
	{
		super(g, args);
		blockN = g.g.cells.size();
		blockSize = grid.A.getN()/blockN;//g.g.cells.get(0).shapeFunctions.size()*g.directions.length;
		blockMatrix = new BlockMatrix(grid.A,blockSize);
		blockDiag = new BlockDiagonalMatrix(grid.A,blockSize);
		System.out.println("inverting");
		blockInverse = blockDiag.inverse();
		System.out.println("inverse");
		downWindOrderings = new ArrayList[grid.directions.length];
		for(int i = 0; i < grid.directions.length; i++)
		{
			downWindOrderings[i] = getDownwindOrdering(grid.directions[i]);
		}
	}

	@Override
	public DoubleTensor smooth(DoubleTensor iterate, DoubleTensor rightHandSide)
	{
		for(int j = 0; j < grid.directions.length; j++)
		{
			ArrayList<Integer> downwindOrdering = downWindOrderings[j];
			for(int k = 0; k < blockN; k++)
			{
				int orderedK = downwindOrdering.get(k);
				DoubleTensor newIterate = getCellVector(orderedK,rightHandSide);
				newIterate = newIterate.sub(IntStream.concat(IntStream.range(0,k),
					IntStream.range(k+1,blockN)).mapToObj(i->
				getCellMatrix(orderedK,downwindOrdering.get(i),blockMatrix).mvmul(getCellVector(downwindOrdering.get(i),
					iterate))).reduce(DoubleTensor::add).orElse(new DoubleTensor(newIterate.size())));

				newIterate = getCellMatrix(orderedK,orderedK,blockInverse).mvmul(newIterate);
				//newIterate.print_formatted();
				setCellVector(orderedK,iterate,newIterate);
			}
		}

		return iterate;
	}
}
