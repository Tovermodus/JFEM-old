import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.sparse.BiCGstab;
import no.uib.cipr.matrix.sparse.GMRES;
import no.uib.cipr.matrix.sparse.IterativeSolverNotConvergedException;
import no.uib.cipr.matrix.sparse.LinkedSparseMatrix;
import org.ojalgo.matrix.store.Primitive64Store;
import org.ojalgo.matrix.store.SparseStore;
import org.ujmp.core.DenseMatrix;
import org.ujmp.core.Matrix;

import java.util.ArrayList;

public class DoubleTensor
{
	double [] denseEntries;
	double [] sparseValues;
	int [] sparseXs;
	int [] sparseYs;
	int sparseEntries;
	int cols;
	int rows;
	boolean isSparse;
	int dimension;
	public DoubleTensor(int size)
	{
		dimension = 1;
		rows = size;
		cols = 1;
		this.isSparse=false;
		denseEntries = new double[size];
	}
	public DoubleTensor(int sizeX, int sizeY, boolean isSparse)
	{
		dimension = 2;
		this.isSparse = isSparse;
		if(!isSparse)
			denseEntries = new double[sizeX*sizeY];
		this.cols = sizeX;
		this.rows = sizeY;
		if(isSparse)
		{
			sparseValues = new double[1];
			sparseYs = new int[1];
			sparseXs = new int[1];
			this.sparseEntries = 1;
		}
	}
	public DoubleTensor(DoubleTensor t)
	{
		cols = t.cols;
		rows = t.rows;
		if(t.isSparse)
		{
			sparseXs = t.sparseXs.clone();
			sparseYs = t.sparseYs.clone();
			sparseValues = t.sparseValues.clone();
		}
		else
			denseEntries = t.denseEntries.clone();
		isSparse = t.isSparse;
		dimension = t.dimension;
	}
	public Matrix getUJMPMatrix()
	{
		Matrix mat = DenseMatrix.Factory.zeros(rows,cols);
		if(!isSparse)
			for(int i = 0; i < getM(); i++)
				for(int j = 0; j < getN(); j++)
					mat.setAsDouble(denseEntries[i*getN()+j], i, j);
		else

			for(int i = 0; i < sparseValues.length; i++)
			{
				mat.setAsDouble(mat.getAsDouble(sparseYs[i],sparseXs[i])+sparseValues[i],
					sparseYs[i],
					sparseXs[i]);
			}
		return mat;
	}
	private SparseStore<Double> getOJALGOMatrix()
	{

		if(dimension != 2|| !isSparse)
			throw new UnsupportedOperationException();
		SparseStore<Double> m = SparseStore.PRIMITIVE64.make(getM(),getN());
		for(int i = 0; i < sparseValues.length; i++)
		{
			m.add(sparseYs[i],sparseXs[i],sparseValues[i]);
		}
		return m;
	}
	private LinkedSparseMatrix getMTJMatrix()
	{
		if(dimension != 2|| !isSparse)
			throw new UnsupportedOperationException();
		LinkedSparseMatrix m = new LinkedSparseMatrix(getM(),getN());
		for(int i = 0; i < sparseValues.length; i++)
		{
			m.add(sparseYs[i],sparseXs[i],sparseValues[i]);
		}
		return m;
	}
	public Matrix getUJMPVector()
	{
		Matrix mat = DenseMatrix.Factory.zeros(rows,1);
		for(int i = 0; i < size(); i++)
			mat.setAsDouble(denseEntries[i], i, 0);
		return mat;
	}
	private Primitive64Store getOJALGOVector()
	{
		if(dimension != 1)
			throw new UnsupportedOperationException();
		Primitive64Store m = Primitive64Store.FACTORY.make(size(),1);
		for(int i = 0; i < size(); i++)
		{
			m.set(i,at(i));
		}
		return m;
	}
	private DenseVector getMTJvector()
	{
		if(dimension != 1)
			throw new UnsupportedOperationException();
		DenseVector m = new DenseVector(size());
		for(int i = 0; i < size(); i++)
		{
			m.set(i,at(i));
		}
		return m;
	}
	public final DoubleTensor fromUJMPVector(Matrix vector)
	{
		DoubleTensor ret = new DoubleTensor((int)vector.getRowCount());
		for(int i = 0; i < ret.size(); i++)
		{
			ret.denseEntries[i] = vector.getAsDouble(i,0);
		}
		return ret;
	}
	private DoubleTensor fromOJALGOVector(Primitive64Store v)
	{
		DoubleTensor m = new DoubleTensor(v.size());
		for(int i = 0; i < v.size(); i++)
		{
			m.set(i,v.get(i));
		}
		return m;
	}
	private DoubleTensor fromMTJvector(DenseVector v)
	{
		DoubleTensor m = new DoubleTensor(v.size());
		for(int i = 0; i < v.size(); i++)
		{
			m.set(i,v.get(i));
		}
		return m;
	}
	public int size()
	{
		return rows* cols;
	}
	public int getM()
	{
		if(dimension != 2)
			throw new UnsupportedOperationException();
		return rows;
	}
	public int getN()
	{
		if(dimension != 2)
			throw new UnsupportedOperationException();
		return cols;
	}
	public double at(int vectorIndex)
	{
		if(dimension != 1)
			throw new UnsupportedOperationException();
		return denseEntries[vectorIndex];
	}
	public double at(int matrixIndex1, int matrixIndex2)
	{
		if(dimension != 2)
			throw new UnsupportedOperationException();
		if(!isSparse)
			return denseEntries[matrixIndex1*getN()+matrixIndex2];
		else
		{
			double ret = 0;
			for(int i = 0; i < sparseEntries; i++)
			{
				if(sparseYs[i] == matrixIndex1)
					if(sparseXs[i] == matrixIndex2)
						ret += sparseValues[i];
			}
			return ret;
		}
	}
	public void set(int vectorIndex, double value)
	{
		if(dimension != 1)
			throw new UnsupportedOperationException();
		denseEntries[vectorIndex] = value;
	}
	private void resizeSparse()
	{
		double [] sparseVals = new double[sparseValues.length*2];
		int [] sparseX  = new int[sparseValues.length*2];
		int [] sparseY =  new int[sparseValues.length*2];
		for(int i = 0; i < sparseEntries; i++)
		{
			sparseVals[i] = sparseValues[i];
			sparseX[i] = sparseXs[i];
			sparseY[i] = sparseYs[i];
		}
		sparseValues = sparseVals;
		sparseYs = sparseY;
		sparseXs = sparseX;
	}

	public void set(int matrixIndex1,int matrixIndex2, double value)
	{
		if(dimension != 2)
			throw new UnsupportedOperationException();
		if(!isSparse)
			denseEntries[matrixIndex1*getN()+matrixIndex2] = value;
		else
		{
			if(sparseEntries == sparseValues.length)
				resizeSparse();
			sparseValues[sparseEntries] = value - at(matrixIndex1,matrixIndex2);
			sparseYs[sparseEntries] = matrixIndex1;
			sparseXs[sparseEntries] = matrixIndex2;
			sparseEntries++;
		}

	}
	public void add(int matrixIndex1,int matrixIndex2, double value)
	{
		if(dimension != 2)
			throw new UnsupportedOperationException();
		if(isSparse)
		{
			if(sparseEntries == sparseValues.length)
				resizeSparse();
			sparseValues[sparseEntries] = value;
			sparseYs[sparseEntries] = matrixIndex1;
			sparseXs[sparseEntries] = matrixIndex2;
			sparseEntries++;
		}
		else
			denseEntries[matrixIndex1*getN()+matrixIndex2]+= value;
	}
	public void add(int vectorIndex, double value)
	{
		if(dimension != 1)
			throw new UnsupportedOperationException();
		denseEntries[vectorIndex]+= value;
	}
	public void print_formatted()
	{
		System.out.println();
		if(dimension == 2)
			for(int i = 0; i < getM(); i++)
			{
				System.out.print("[");
				for(int j = 0; j < getN(); j++)
				{
					if (Math.abs(at(i, j)) <= 1e-14)
						System.out.print("<>........  ");
					else
					{
						if(at(i,j)>=0)
							System.out.print("+");
						System.out.print(String.format("%6.3e", at(i, j)) + "  ");
					}
				}
				System.out.println("]");
			}
		if(dimension == 1)
		{
			System.out.print("[");
			for (int j = 0; j < size(); j++)
			{
				if (Math.abs(at(j)) <= 1e-14)
					System.out.print("<>........  ");
				else
				{
					if(at(j)>=0)
						System.out.print("+");
					System.out.print(String.format("%6.3e", at(j)) + "  ");
				}
			}
			System.out.println("]");
		}
		System.out.println();
	}

	public void print_formatted(String s)
	{
		System.out.println("   "+s);
		if(dimension == 2)
			for(int i = 0; i < getM(); i++)
			{
				System.out.print("[");
				for(int j = 0; j < getN(); j++)
				{
					if (Math.abs(at(i, j)) <= 1e-14)
						System.out.print("<>........  ");
					else
					{
						if(at(i,j)>=0)
							System.out.print("+");
						System.out.print(String.format("%6.3e", at(i, j)) + "  ");
					}
				}
				System.out.println("]");
			}
		if(dimension == 1)
		{
			System.out.print("[");
			for (int j = 0; j < size(); j++)
			{
				if (Math.abs(at(j)) <= 1e-14)
					System.out.print("<>........  ");
				else
				{
					if(at(j)>=0)
						System.out.print("+");
					System.out.print(String.format("%6.3e", at(j)) + "  ");
				}
			}
			System.out.println("]");
		}
		System.out.println();
	}
	public void print_log()
	{
		if(dimension != 2)
			return;
		for(int i = 0; i < getM(); i++)
		{
			for(int j = 0; j < getN(); j++)
			{
				if (Math.abs(at(i, j)) <= 1e-14)
					System.out.print("   \t");
				else
					System.out.print((int) Math.log10(Math.abs(at(j, i))) + "\t");
			}
			System.out.println();
		}
	}
	public static DoubleTensor vectorFromValues(double... values)
	{

		DoubleTensor ret = new DoubleTensor(values.length);
		for(int i = 0; i < values.length; i++)
			ret.set(i,values[i]);
		return ret;
	}
	public static DoubleTensor squareMatrixFromValues(double... values)
	{
		if(Math.pow((int)Math.sqrt(values.length),2) != values.length)
			throw new UnsupportedOperationException("nope");
		DoubleTensor ret = new DoubleTensor((int)Math.sqrt(values.length),(int)Math.sqrt(values.length),false);
		ret.denseEntries = values;
		return ret;
	}
	public double x()
	{
		if(dimension == 1)
			return at(0);
		else throw new UnsupportedOperationException("nope");
	}
	public double y()
	{
		if(dimension == 1)
			return at(1);
		else throw new UnsupportedOperationException("nope");
	}
	public double inner(DoubleTensor other)
	{
		double ret = 0;
		if(this.size() == other.size())
		{
			for(int i = 0; i < this.size(); i++)
			{
				ret += denseEntries[i]*other.denseEntries[i];
			}
			return ret;
		}
		throw new UnsupportedOperationException();
	}
	public DoubleTensor add(DoubleTensor other)
	{
		if(this.dimension == other.dimension && this.size() == other.size())
		{
			DoubleTensor ret =  new DoubleTensor(this);
			ret.dimension = this.dimension;
			ret.isSparse = this.isSparse && other.isSparse;
			if(!this.isSparse && !other.isSparse)
			{
				for(int i = 0; i < denseEntries.length; i++)
					ret.denseEntries[i] = this.denseEntries[i]+other.denseEntries[i];
			}
			if(!this.isSparse && other.isSparse)
			{
				ret.denseEntries = this.denseEntries;
				for(int i = 0; i < other.sparseEntries; i++)
				{
					ret.denseEntries[other.sparseYs[i]*getN() + other.sparseXs[i]] += other.sparseValues[i];
				}
			}
			if(this.isSparse && !other.isSparse)
			{
				ret.denseEntries = other.denseEntries;
				for(int i = 0; i < this.sparseEntries; i++)
				{
					ret.denseEntries[this.sparseYs[i]*getN() + this.sparseXs[i]] += this.sparseValues[i];
				}
			}
			if(this.isSparse && other.isSparse)
			{
				ret.sparseEntries = this.sparseEntries + other.sparseEntries;

				for(int i = 0; i < this.sparseEntries; i++)
				{
					ret.sparseXs[i] = this.sparseXs[i];
					ret.sparseYs[i] = this.sparseYs[i];
					ret.sparseValues[i] = this.sparseValues[i];
				}
				for(int i = this.sparseEntries; i < other.sparseEntries+this.sparseEntries; i++)
				{
					ret.sparseXs[i] = other.sparseXs[i];
					ret.sparseYs[i] = other.sparseYs[i];
					ret.sparseValues[i] = other.sparseValues[i];
				}
			}
			return ret;
		}
		throw new UnsupportedOperationException();
	}
	public DoubleTensor mul(double scalar)
	{
		DoubleTensor ret = new DoubleTensor(this);
		if(!isSparse)
		{
			for(int i = 0; i < denseEntries.length; i++)
			{
				ret.denseEntries[i] = denseEntries[i]*scalar;
			}
		}
		else
		{
			for(int i = 0; i < sparseEntries; i++)
			{
				ret.sparseValues[i] = sparseValues[i]*scalar;
			}
		}
		ret.dimension = dimension;
		ret.isSparse = isSparse;
		return ret;
	}
	public DoubleTensor mvmul(DoubleTensor vector)
	{
		if(vector.dimension != 1 || this.dimension != 2)
			throw new UnsupportedOperationException();
		DoubleTensor ret = new DoubleTensor(vector.size());
		ret.isSparse = false;
		ret.dimension = 1;
		if(!isSparse)
		{
			for(int i = 0; i < getM(); i++)
			{
				for(int j = 0; j < getN(); j++)
					ret.denseEntries[i] += denseEntries[i*getN()+j]*vector.denseEntries[j];
			}
		}
		else
		{
			for(int i = 0; i < sparseEntries; i++)
			{
				ret.denseEntries[sparseYs[i]]+=sparseValues[i]*vector.denseEntries[sparseXs[i]];
			}
		}
		return ret;
	}

	public DoubleTensor mmmul(DoubleTensor other)
	{
		if(other.dimension != 2 || this.dimension != 2)
			throw new UnsupportedOperationException();
		DoubleTensor ret = new DoubleTensor(this.getM(),other.getN(),this.isSparse && other.isSparse);
		if(!this.isSparse && !other.isSparse)
		{
			for(int i = 0; i < getM(); i++)
				for(int j = 0; j < other.getN(); j++)
					for(int k = 0; k < getN(); k++)
						ret.denseEntries[i*other.getN()+j] +=
							denseEntries[i*getN()+k]*other.denseEntries[k*other.getN()+j];
		}
		if(!this.isSparse && other.isSparse)
		{
			for(int i = 0; i < other.sparseEntries; i++)
				for(int j = 0; j < getM(); j++)
				{
					ret.denseEntries[j * ret.getN() + other.sparseXs[i]] +=
						denseEntries[j * getN() + other.sparseYs[i]] * other.sparseValues[i];
				}
		}
		if(this.isSparse && !other.isSparse)
		{
			for(int i = 0; i < this.sparseEntries; i++)
				for(int j = 0; j < other.getN(); j++)
					ret.denseEntries[sparseYs[i]*ret.getN() + j] +=
						sparseValues[i] * other.denseEntries[sparseXs[i]*other.getN()+j];
		}
		if(this.isSparse && other.isSparse)
		{
			for(int i = 0; i < this.sparseEntries; i++)
				for(int j = 0; j < other.sparseEntries; j++)
				{
					if(this.sparseXs[i] == other.sparseYs[j])
						ret.add(this.sparseYs[i], other.sparseXs[j],
							this.sparseValues[i]*other.sparseValues[j]);
				}
		}
		return ret;
	}
	public DoubleTensor transpose()
	{
		if(this.dimension != 2)
			throw new UnsupportedOperationException();
		DoubleTensor ret = new DoubleTensor(this);
		ret.dimension = 2;
		if(isSparse)
		{
			ret.sparseXs = sparseYs;
			ret.sparseYs = sparseXs;
			ret.sparseValues = sparseValues;
			ret.sparseEntries = sparseEntries;
		}
		else
		{
			for(int i = 0; i < getM(); i++)
				for(int j = 0; j < getN(); j++)
				{
					ret.denseEntries[j*getM()+i] = denseEntries[i*getN()+j];
				}
		}
		ret.rows = cols;
		ret.cols = rows;
		ret.isSparse = isSparse;
		return ret;
	}
	public double vectorNorm()
	{
		if(this.dimension != 1)
			throw new UnsupportedOperationException();
		return Math.sqrt(this.inner(this));
	}
	public DoubleTensor sub(DoubleTensor other)
	{
		if(this.dimension != other.dimension && this.size() != other.size())
		{
			throw new UnsupportedOperationException();
		}
		return this.add(other.mul(-1.));
	}
	public DoubleTensor solveSymm(DoubleTensor rightHandSide)
	{
		if(rightHandSide.dimension!=1)
			throw new UnsupportedOperationException();
		DoubleTensor ret;
		ret = fromUJMPVector(this.getUJMPMatrix().solveSymm(rightHandSide.getUJMPVector()));
		return ret;
	}
	public DoubleTensor solve(DoubleTensor rightHandSide)
	{
		if(rightHandSide.dimension!=1)
			throw new UnsupportedOperationException();
		DoubleTensor ret;
		ret = fromUJMPVector(this.getUJMPMatrix().solve(rightHandSide.getUJMPVector()));
		ret.print_formatted();
		return ret;
	}
	public DoubleTensor solveCG(DoubleTensor rightHandSide, double residualTolerance)
	{
		//System.out.println(this.mat.isSparse());

		if(getN() != getN() || dimension != 2 || rightHandSide.dimension != 1)
			throw new UnsupportedOperationException();
		DoubleTensor z;
		DoubleTensor newResiduum;
		DoubleTensor iterate = new DoubleTensor(getM());
		DoubleTensor residuum = rightHandSide.sub(this.mvmul(iterate));
		DoubleTensor defect = new DoubleTensor(residuum);
		double α = 0;
		double β = 0;
		for(int iter = 0; iter < getM() && residuum.vectorNorm() > residualTolerance; iter++)
		{
			z = this.mvmul(defect);
			α = residuum.inner(residuum)/defect.inner(z);
			iterate = iterate.add(defect.mul(α));
			newResiduum = residuum.sub(z.mul(α));
			β = newResiduum.inner(newResiduum)/residuum.inner(residuum);
			defect = newResiduum.add(defect.mul(β));
			residuum = newResiduum;
			if(iter%10 == 0)
				System.out.println("iter "+iter+" " +residuum.vectorNorm());
		}
		return iterate;
	}
	public DoubleTensor solveGMRES(DoubleTensor rightHandSide)
	{
		DenseVector rhs = rightHandSide.getMTJvector();
		DenseVector v = new DenseVector(getM());
		DenseVector x = new DenseVector(getM());
		LinkedSparseMatrix m = getMTJMatrix();
		GMRES gm = new GMRES(v);
		try
		{
			x = (DenseVector) gm.solve(m, rhs, x);
		} catch (IterativeSolverNotConvergedException e)
		{
			e.printStackTrace();
		}
		System.out.println(x.size());
		return fromMTJvector(x);
	}
	public DoubleTensor solveBiCGStab(DoubleTensor rightHandSide, double residualTolerance)
	{
		if(getN() != getN() || dimension != 2 || rightHandSide.dimension != 1)
			throw new UnsupportedOperationException();
		DoubleTensor v;
		DoubleTensor s;
		DoubleTensor t;
		DoubleTensor iterate = new DoubleTensor(getM());
		DoubleTensor startResiduum = rightHandSide.sub(this.mvmul(iterate));
		DoubleTensor residuum = new DoubleTensor(startResiduum);
		DoubleTensor p = new DoubleTensor(residuum);
		double alpha;
		double omega;
		double rho = residuum.inner(residuum);
		double rhoLast;
		double beta;
		for(int iter = 0; iter < getM() && residuum.vectorNorm() > residualTolerance; iter++)
		{
			v = this.mvmul(p);
			alpha = rho/v.inner(startResiduum);
			s = residuum.sub(v.mul(alpha));
			t = this.mvmul(s);
			omega = t.inner(s)/t.inner(t);
			iterate = iterate.add(p.mul(alpha)).add(s.mul(omega));
			residuum = s.sub(t.mul(omega));
			rhoLast = rho;
			rho = residuum.inner(startResiduum);
			beta = alpha/omega*rho/rhoLast;
			p = residuum.add(p.mul(beta)).sub(v.mul(omega*beta));
			if(iter%10 == 0)
				System.out.println("iter "+iter+" "+residuum.vectorNorm());
		}
		return iterate;


	}
}
