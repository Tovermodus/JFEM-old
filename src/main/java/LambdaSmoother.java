import com.google.common.collect.Iterables;

import java.util.stream.IntStream;

public class LambdaSmoother extends AngularSmoother
{
	double isoScatter;
	BlockDiagonalMatrix blockDiag;
	BlockDiagonalMatrix blockInverse;

	public LambdaSmoother(AngularGrid g, String[] args)
	{
		super(g, args);
		int blockSize;
		blockSize = g.g.shapeFunctions.size();
		System.out.println("create blockdiag");
		blockDiag = new BlockDiagonalMatrix(g.Areshuffled,blockSize);
		System.out.println("inverting");
		//blockInverse = blockDiag.inverse();
		System.out.println("inverse");
	}

	@Override
	public DoubleTensor smooth(DoubleTensor iterate, DoubleTensor rightHandSide)
	{
		DoubleTensor res = rightHandSide.sub(grid.A.mvmul(iterate));
		iterate = iterate.add(grid.shuffle.transpose().mvmul(blockDiag.solve(grid.shuffle.mvmul(res))));

		//DoubleTensor res = rightHandSide.sub(grid.A.mvmul(iterate));
		//iterate = iterate.add(grid.TransAbs.solve(res));
		return iterate;
	}
}
