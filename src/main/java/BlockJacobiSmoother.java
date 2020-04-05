public class BlockJacobiSmoother extends AngularSmoother
{
	BlockDiagonalMatrix blockDiag;
	BlockDiagonalMatrix blockInverse;
	double omega;
	public BlockJacobiSmoother(AngularGrid g, String[] args)
	{
		super(g, args);
		int blockSize = g.g.cells.get(0).shapeFunctions.size()*g.directions.length;
		blockDiag = new BlockDiagonalMatrix(grid.A,blockSize);
		System.out.println("inverting");
		blockInverse = blockDiag.inverse();
		System.out.println("inverse");
	}

	@Override
	public DoubleTensor smooth(DoubleTensor iterate, DoubleTensor rightHandSide)
	{
		DoubleTensor res = rightHandSide.sub(grid.A.mvmul(iterate));
		iterate = iterate.add(blockInverse.mvmul(res));

		return iterate;
	}
}
