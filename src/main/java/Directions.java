public enum Directions
{
	TGLC1(new DoubleTensor[]{DoubleTensor.vectorFromValues(0.577,0.577),
		DoubleTensor.vectorFromValues(-0.577,0.577),
		DoubleTensor.vectorFromValues(0.577,-0.577),
		DoubleTensor.vectorFromValues(-0.577,-0.577)},
		new double[]{0.25,0.25,0.25,0.25}),
	K(new DoubleTensor[]{DoubleTensor.vectorFromValues(1,0),
		DoubleTensor.vectorFromValues(-1,0),
		DoubleTensor.vectorFromValues(0,1),
		DoubleTensor.vectorFromValues(0,-1)},
		new double[]{0.25,0.25,0.25,0.25}),
	KS(new DoubleTensor[]{DoubleTensor.vectorFromValues(1/Math.sqrt(2),1/Math.sqrt(2)),
		DoubleTensor.vectorFromValues(-1/Math.sqrt(2),-1/Math.sqrt(2))},
		new double[]{0.5,0.5}),
	TGLC2(new DoubleTensor[]{DoubleTensor.vectorFromValues(0.869,0.360),
		DoubleTensor.vectorFromValues(0.360,0.869),
		DoubleTensor.vectorFromValues(0.359,0.359),
		DoubleTensor.vectorFromValues(-0.869,0.360),
		DoubleTensor.vectorFromValues(-0.360,0.869),
		DoubleTensor.vectorFromValues(-0.359,0.359),
		DoubleTensor.vectorFromValues(0.869,-0.360),
		DoubleTensor.vectorFromValues(0.360,-0.869),
		DoubleTensor.vectorFromValues(0.359,-0.359),
		DoubleTensor.vectorFromValues(-0.869,-0.360),
		DoubleTensor.vectorFromValues(-0.360,-0.869),
		DoubleTensor.vectorFromValues(-0.359,-0.359)},new double[]{0.0816,0.0816,0.087,0.0816,0.0816,0.087,
		0.0816,0.0816,0.087,0.0816,0.0816,0.087});
	private final DoubleTensor[] directions;
	private final double[] directionWeights;
	public DoubleTensor[] getDirections()
	{
		return directions;
	}
	public double[] getDirectionWeights()
	{
		return directionWeights;
	}

	Directions(DoubleTensor[] directions, double[] directionWeights)
	{
		this.directions = directions;
		this.directionWeights = directionWeights;
	}
}
