/*
 * Bayesian Composite Models
 * 
 * (c) 2015 Collaborative Drug Discovery, Inc.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */

package com.cdd.bayes;

import java.io.*;
import java.lang.*;
import java.util.*;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.*;
import org.openscience.cdk.fingerprint.model.*;
import org.openscience.cdk.interfaces.IAtomContainer;

/*
 * Composite models: the main entrypoint class for taking a collection of molecules & data and partitioning into some
 * number of segments, each of which represents a range of continuous values from the input. Once the composite model
 * has been created, prospective molecules can be submitted, and the relative likelihood of the molecule's property
 * being in one of these binned ranges is made available.
 * 
 * The input assumes that property values are "energy-like" rather than exponential, which in practical terms means
 * that bioactivity concentrations should be converted to the log scale (-log10) before use.
 */
public class CompositeModel
{
	protected int minBins = 3, maxBins = 8; // reasonable defaults
	protected double[] segments = null; // default of null means they should be calculated
	protected double minVal = Double.NaN, maxVal = Double.NaN; // lowest & highest values in training set

	public static final class Entry
	{
		public IAtomContainer mol = null;
		public double val = Double.NaN;
		public int[] fp = null;
	}

	protected List<Entry> entries = new ArrayList<Entry>();

	protected Bayesian[] models = null; // the payload: one model per bin is delivered
	protected int[][] matrix = null; // validation matrix [want][got]: diagonal entries are hits, off-diagonals are miss-by-distance

	// constants used for internal workings; may need to tweak these to get optimal results
	private final int CLUSTER_SUBSIZE = 100; // largest size of subset used for pre-clustering to estimate ROCs (high throughput)
	private final int MAX_CANDIDATES = 50; // number of candidate segments to consider (in order of best first)
	private final float MIN_ROC_SPLIT = 0.55f; // when best ROC for splitting a segment drops below this value, stop
	private final float MIN_BIN_FRACTION = 0.05f; // creating a bin with less than this portion of entries is disallowed

	// ------------ public methods ------------

	// create an empty model: it needs to be filled up with training data then calculated before it can be used for prediction
	public CompositeModel()
	{
	}

	// create a model using data that was previously built; the boundary parameter includes the segmentation breaks,
	// as well as the minimum/maximum values (see getBoundaries()); the models need to have been created with a previous
	// instance (see getModels()); note that the models list is shallow-copied
	public CompositeModel(double[] boundary, Bayesian[] models)
	{
		final int nbins = models.length;

		if (boundary.length < 4 || boundary.length != nbins + 1) 
			throw new ModelException("Invalid parameters: boundary length=" + boundary.length + ", model length=" + nbins);

		minVal = boundary[0];
		maxVal = boundary[boundary.length - 1];
		segments = new double[boundary.length - 2];
		for (int n = 0; n < segments.length; n++) segments[n] = boundary[n + 1];

		for (Bayesian m : models) if (m.getClassType() != CircularFingerprinter.CLASS_ECFP6) 
			throw new ModelException("Only ECFP6 fingerprints are permitted for composite models.");

		this.models = Arrays.copyOf(models, nbins);
	}

	// adds a single molecule & value to the collection of contents that will be operated upon
	public void addEntry(IAtomContainer mol, double val)
	{
		Entry e = new Entry();
		e.mol = mol;
		e.val = val;
		entries.add(e);
	}
	
	// add an already-instantiated entry to the list
	public void addEntry(Entry e)
	{
		fillFingerprints(e);
		entries.add(e);
	}
	
	// makes sure the fingerprint field is defined
	public void fillFingerprints(Entry e)
	{
		if (e.fp != null) return;

		CircularFingerprinter circ = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP6);
		try {circ.calculate(e.mol);}
		catch (CDKException ex) {throw new ModelException(ex);}
		Set<Integer> fplist = new TreeSet<Integer>();
		for (int n = circ.getFPCount() - 1; n >= 0; n--) fplist.add(circ.getFP(n).hashCode);

		e.fp = new int[fplist.size()];
		int p = 0;
		for (int h : fplist) e.fp[p++] = h;
	}

	// access to user-provided molecule/value/fingerprint content
	public int numEntries()
	{
		return entries.size();
	}

	public Entry getEntry(int N)
	{
		return entries.get(N);
	}

	// control over number of bins; note that it is valid for min & max to be the same (reduces degrees of freedom); minimum
	// number of bins is 3, since any less than that defeats the purpose of using the composite model in the first place
	public int getMinBins()
	{
		return minBins;
	}

	public int getMaxBins()
	{
		return maxBins;
	}

	public void setMinBins(int nbins)
	{
		minBins = Math.max(3, nbins);
	}

	public void setMaxBins(int nbins)
	{
		maxBins = Math.min(20, nbins);
	}

	public void setNumBins(int nbins)
	{
		setMinBins(nbins);
		setMaxBins(nbins);
	}

	// segments are the cutpoints for separating the bins, with the idea number being #bins-1; the user may provide any number of
	// them, some or all of which may be used; the cut points will be calculated automatically if not provided
	public double[] getSegments()
	{
		return segments;
	}

	public void setSegments(double[] seg)
	{
		segments = seg;
	}

	// perform the calculation: assuming that entries have been provided, and any other preparation has been done, proceed to
	// create all of the bins, sub-models, and calibrations
	public void calculate() throws CDKException
	{
		final int num = entries.size();
		if (num == 0) throw new ModelException("No entries provided.");
		if (num < minBins) throw new ModelException("Min bins=" + minBins + " and # entries=" + num + ": this isn't going to work.");
		if (segments != null && segments.length >= num - 1) 
			throw new ModelException("Provided " + num + " entries and " + segments.length + " segments: this isn't going to work.");

		if (segments == null) determineSegments();

		// record min/max
		minVal = Double.POSITIVE_INFINITY;
		maxVal = Double.NEGATIVE_INFINITY;
		for (Entry e : entries)
		{
			minVal = Math.min(minVal, e.val);
			maxVal = Math.max(maxVal, e.val);
		}

		// prepare bin assignments
		int[][] bins = assignBins(segments);
		final int nbins = bins.length;
		int[] binidx = new int[num];
		for (int n = 0; n < nbins; n++) for (int b : bins[n]) binidx[b] = n;

		// generate a calibrated model for each bin
		models = new Bayesian[nbins];
		for (int n = 0; n < nbins; n++)
		{
			models[n] = new Bayesian(CircularFingerprinter.CLASS_ECFP6,0);
			for (int i = 0; i < nbins; i++) for (int b : bins[i])
				models[n].addMolecule(entries.get(b).mol, i == n);
			models[n].build();
			models[n].validateFiveFold();
			//Main.writeln("  bin="+n+" roc="+models[n].getROCAUC());
		}

		// validation matrix: mapping is [want][got]
		matrix = new int[nbins][];
		for (int n = 0; n < nbins; n++) matrix[n] = new int[nbins];

		for (int n = 0; n < num; n++)
		{
			int best = -1;
			double highest = Double.NEGATIVE_INFINITY;
			for (int i = 0; i < nbins; i++)
			{
				IAtomContainer mol = entries.get(n).mol;
				double v = models[i].scalePredictor(models[i].predict(mol));
				if (v > highest)
				{
					best = i;
					highest = v;
				}
			}
			matrix[binidx[n]][best]++;
		}

		//for (int n=0;n<nbins;n++) Main.writeln(Arrays.toString(matrix[n]));
	}

	// performs an automated determination of viable segments - the cutpoints for binning - based on the entries that have
	// been provided; normally this is called by the calculate() method above, but it can be called separately preemptively,
	// if for some reason the caller wants to inspect and modify
	public void determineSegments() throws CDKException
	{
		final int num = entries.size();
		if (num == 0) throw new ModelException("No entries provided.");

		// obtain a reasonable subset: this should be small enough that building a model for every possible permutation is
		// not a rate limiting performance issue
		int[] subset;
		if (num > CLUSTER_SUBSIZE)
		{
			GreedyLinearCluster glc = new GreedyLinearCluster(entries,CLUSTER_SUBSIZE);
			subset = glc.calculate();
		}
		else
		{
			subset = new int[num];
			for (int n = 0; n < num; n++) subset[n] = n;
		}

		// enumerate all valid interstitial cutpoints; each of these is technically fair game for proposal as one of the
		// segment boundaries
		final int sz = subset.length;
		float[] values = new float[sz];
		for (int n = 0; n < sz; n++) values[n] = (float) entries.get(subset[n]).val;
		Arrays.sort(values);
		List<Double> cuts = new ArrayList<Double>();
		for (int n = 0; n < sz - 1; n++)
			if (values[n] != values[n + 1]) cuts.add(0.5 * (values[n] + values[n + 1]));
		int ncuts = cuts.size();

		if (ncuts < minBins) throw new ModelException("Unable to find reasonable number of cut points.");

		// compute a ROC score for each of the putative cutpoints
		float[] integrals = new float[ncuts];
		for (int n = 0; n < ncuts; n++) integrals[n] = sampleBayesianROC(subset, cuts.get(n));
		float imin = integrals[0], imax = imin;
		for (int n = 1; n < ncuts; n++)
		{
			imin = Math.min(imin, integrals[n]);
			imax = Math.max(imax, integrals[n]);
		}
		float iscale = 1 / (imax - imin);
		for (int n = 0; n < ncuts; n++) integrals[n] = (integrals[n] - imin) * iscale;

		// plot all points as Gaussians, to be able to measure their height and gradient
		final int npt = 1000;
		float lowV = values[0], highV = lowV;
		for (int n = 1; n < values.length; n++)
		{
			lowV = Math.min(lowV, values[n]);
			highV = Math.max(highV, values[n]);
		}
		lowV -= 1;
		highV += 1;
		final float gauss = 0.1f / (highV - lowV) * npt;
		float[] height = new float[npt];
		for (float v : values) plotGaussian(height, (v - lowV) / (highV - lowV) * npt, gauss);

		// figure out the 2nd derivative: higher values are for local minima, which are good
		float[] deriv1 = new float[npt], deriv2 = new float[npt];
		for (int n = 1; n < npt - 1; n++) deriv1[n] = height[n + 1] - height[n - 1]; // (no need to calibrate)
		for (int n = 2; n < npt - 2; n++) deriv2[n] = deriv1[n + 1] - deriv1[n - 1]; // (ditto)
		float dmin = deriv2[0], dmax = dmin;
		for (int n = 1; n < npt; n++)
		{
			dmin = Math.min(dmin, deriv2[n]);
			dmax = Math.max(dmax, deriv2[n]);
		}
		float dscale = 1 / (dmax - dmin);
		for (int n = 0; n < npt; n++) deriv2[n] = (deriv2[n] - dmin) * dscale;

		// for each of the cutpoints, add up the area above & below
		float[] ratio = new float[ncuts];
		int minBinSize = (int) Math.ceil(MIN_BIN_FRACTION * entries.size());
		for (int n = 0; n < ncuts; n++)
		{
			int above = 0, below = 0;
			final float cutval = cuts.get(n).floatValue();
			for (Entry e : entries)
			{
				if (e.val >= cutval) above++;
				else below++;
			}

			if (above < minBinSize || below < minBinSize) // rare, but it happens
			{
				cuts.remove(n);
				ncuts--;
				n--;
				continue;
			}

			ratio[n] = Math.max((above + 1.0f) / (below + 1.0f), (below + 1.0f) / (above + 1.0f));
		}
		if (ncuts < minBins) throw new ModelException("Unable to find reasonable number of cut points.");

		float[] sortedRatio = Arrays.copyOf(ratio, ncuts);
		Arrays.sort(sortedRatio);
		float ratioMax = sortedRatio[(int) (0.9f * ncuts)], rscale = 1 / ratioMax; // scale so that 1 --> 90% by popularity
		for (int n = 0; n < ncuts; n++) ratio[n] *= rscale;

		// now calculate a desirability for each of the cutpoints; lower is better
		final float[] desire = new float[ncuts];
		for (int n = 0; n < ncuts; n++)
		{
			// start with ROC integrals; prescaled 0=worst, 1=best
			desire[n] += 1 - integrals[n];

			// add in scaled 2nd derivative, where close to 0 corresponds to a saddle or local minimum
			float px = (cuts.get(n).floatValue() - lowV) / (highV - lowV) * npt;
			desire[n] += 1 - interpolate(deriv2, px);

			// add in ratio: perfect balance is 1, way off the end tends toward N
			desire[n] += ratio[n];
		}

		// for debugging purposes only
		/*Main.writeln("Cutoffs ("+ncuts+")");
		for (int n=0;n<ncuts;n++) Main.writeln("  value="+cuts.get(n).floatValue()+" desirability="+desire[n]);*/

		// pick the best ones: sort by desirability, and pick the best one as the primary segmentation boundary; the next
		// few can be considered as possibilities for the subsequent iterative refinement
		Integer[] sorted = new Integer[ncuts];
		for (int n = 0; n < ncuts; n++) sorted[n] = n;
		Arrays.sort(sorted, new Comparator<Integer>()
		{
			public int compare(Integer i1, Integer i2)
			{
				final float v1 = desire[i1], v2 = desire[i2];
				if (v1 < v2) return -1;
				else if (v1 > v2) return 1;
				else return 0;
			}
		});
		segments = new double[]{cuts.get(sorted[0])};
		List<Double> candidates = new ArrayList<Double>();
		for (int n = 1; n < ncuts && n < MAX_CANDIDATES; n++) candidates.add(cuts.get(sorted[n]));

		iterativelyAddSegments(candidates);
	}

	// obtain information about models and validation
	public int numBins()
	{
		return segments == null ? 0 : segments.length + 1;
	}

	public int[][] getAssignedBins()
	{
		return assignBins(segments);
	}

	public Bayesian getModel(int N)
	{
		return models[N];
	}

	public Bayesian[] getModels()
	{
		return models;
	}

	public int[][] getValidationMatrix()
	{
		return matrix;
	}

	// the "boundaries" is an array that includes the segments, capped by the minimum & maximum values; its
	// size is #bins + 1; storing the min/max value affects not at all the model creation or application, but
	// it can be useful for interpreting the assignments
	public double getMinVal()
	{
		return minVal;
	}

	public double getMaxVal()
	{
		return maxVal;
	}

	public double[] getBoundaries()
	{
		double[] bound = new double[segments.length + 2];
		bound[0] = minVal;
		for (int n = 0; n < segments.length; n++) bound[n + 1] = segments[n];
		bound[segments.length + 1] = maxVal;
		return bound;
	}

	// using the model to make new predictions; the result is an array with calibrated prediction scores for each of the
	// available bins, whereby most values should be in the range of (0..1); the highest value can be considered to be the
	// winner, but other bins with comparable scores might be contenders
	public float[] predictBins(IAtomContainer mol) throws CDKException
	{
		final int nbins = models.length;
		float[] pred = new float[nbins];
		for (int n = 0; n < nbins; n++)
		{
			pred[n] = (float) models[n].scalePredictor(models[n].predict(mol));
		}
		return pred;
	}

	// ------------ private methods ------------

	// for the given subset (by index) and threshold for activity, build a Bayesian model and return its ROC integral
	private float sampleBayesianROC(int[] subset, double threshold) throws CDKException
	{
		Bayesian bayes = new Bayesian(CircularFingerprinter.CLASS_ECFP6,0);
		for (int i : subset)
		{
			Entry e = entries.get(i);
			bayes.addMolecule(e.mol, e.val >= threshold);
		}
		bayes.build();
		bayes.validateLeaveOneOut();
		return (float) bayes.getROCAUC();
	}

	// similar to above, but takes two pre-formed partitions as arbitrary true/false
	private float sampleBayesianROC(List<Entry> ptn1, List<Entry> ptn2) throws CDKException
	{
		Bayesian bayes = new Bayesian(CircularFingerprinter.CLASS_ECFP6,0);
		for (Entry e : ptn1) bayes.addMolecule(e.mol, false);
		for (Entry e : ptn2) bayes.addMolecule(e.mol, true);
		bayes.build();
		bayes.validateLeaveOneOut();
		return (float) bayes.getROCAUC();
	}

	// increments the Y values based on a gaussian distribution at X, which is chosen to have an integral of 1
	private void plotGaussian(final float[] y, float x, float sigma)
	{
		final float k = 1.0f / (float) Math.sqrt(Math.PI * 2 * sigma * sigma); // scales the integral to 1
		final float a = 1.0f / (2 * sigma * sigma);
		float total = 0;
		for (int n = 0; n < y.length; n++)
		{
			final float d = n - x;
			final float v = k * (float) Math.exp(-d * d * a);
			total += v;
			y[n] += v;
		}

		if (total < 0.99f)
		{
			float extra = 1 - total;
			if (2 * x < y.length) y[0] += extra;
			else y[y.length - 1] += extra;
		}
	}

	// picks an implied y-value for an off-integral point
	private float interpolate(float[] y, float x)
	{
		if (x < 0) return y[0];
		if (x >= y.length - 1) return y[y.length - 1];
		int ix = (int) Math.floor(x);
		float rx = x - ix;
		return y[ix] * (1 - rx) + y[ix + 1] * rx;
	}

	// given that there is at least one segment already defined, contemplate adding a number of additional segments, and
	// evaluate the effects of making these into bins
	private void iterativelyAddSegments(List<Double> candidates) throws CDKException
	{
		if (candidates.size() == 0 || segments.length >= maxBins - 1) return;

		int minBinSize = (int) Math.ceil(MIN_BIN_FRACTION * entries.size());

		int bestCandidate = -1;
		double[] bestNewSeg = null;
		float bestROC = 0;

		//Main.writeln("Iteratively considering:"+candidates);

		for (int n = 0; n < candidates.size(); n++)
		{
			double cseg = candidates.get(n);
			int idx = segments.length;
			double[] newseg = Arrays.copyOf(segments, idx + 1);
			for (; idx > 0 && newseg[idx - 1] > cseg; idx--) newseg[idx] = newseg[idx - 1];
			newseg[idx] = cseg;

			// partition the bins with the new segment, and check to see if any of them cause a too-small partition
			int[][] bins = assignBins(newseg);
			/*Main.write("  n="+n+" cseg="+cseg+" newseg="+Arrays.toString(newseg)+" bins=");
			for (int i=0;i<bins.length;i++) Main.write((i>0 ? "," : "")+bins[i].length);
			Main.writeln();*/
			boolean anyTooSmall = false;
			for (int[] b : bins) if (b.length < minBinSize)
			{
				anyTooSmall = true;
				break;
			}
			if (anyTooSmall)
			{
				candidates.remove(n);
				n--;
				continue;
			}

			// look at the two new bins that were created by adding the segment, and see how well they separate
			// from each other, by creating a tentative model
			List<Entry> ptn1 = new ArrayList<Entry>(), ptn2 = new ArrayList<Entry>();
			for (int i : bins[idx]) ptn1.add(entries.get(i));
			for (int i : bins[idx + 1]) ptn2.add(entries.get(i));
			if (ptn1.size() > CLUSTER_SUBSIZE)
			{
				GreedyLinearCluster glc = new GreedyLinearCluster(ptn1,CLUSTER_SUBSIZE);
				List<Entry> ptn = new ArrayList<Entry>();
				for (int i : glc.calculate()) ptn.add(ptn1.get(i));
				ptn1 = ptn;
			}
			if (ptn2.size() > CLUSTER_SUBSIZE)
			{
				GreedyLinearCluster glc = new GreedyLinearCluster(ptn2,CLUSTER_SUBSIZE);
				List<Entry> ptn = new ArrayList<Entry>();
				for (int i : glc.calculate()) ptn.add(ptn2.get(i));
				ptn1 = ptn;
			}
			float roc = sampleBayesianROC(ptn1, ptn2);

			if (bestCandidate < 0 || roc > bestROC)
			{
				bestCandidate = n;
				bestNewSeg = newseg;
				bestROC = roc;
			}
		}

		if (bestCandidate < 0) return;
		if (segments.length > 1 && bestROC < MIN_ROC_SPLIT) return;
		segments = bestNewSeg;
		candidates.remove(bestCandidate);
		iterativelyAddSegments(candidates);
	}

	// given a set of putative segment boundaries, makes a list of bins and the entries that fall into them
	private int[][] assignBins(double[] seg)
	{
		final int nbins = seg.length + 1, nent = entries.size();
		int[] binidx = new int[nent], binsz = new int[nbins];
		for (int n = 0; n < nent; n++)
		{
			final double v = entries.get(n).val;
			for (int i = 0; i < seg.length; i++) if (v >= seg[i]) binidx[n] = i + 1;
			binsz[binidx[n]]++;
		}

		int[][] bins = new int[nbins][];
		for (int n = 0; n < nbins; n++)
		{
			bins[n] = new int[binsz[n]];
			binsz[n] = 0;
		}
		for (int n = 0; n < nent; n++)
		{
			final int b = binidx[n];
			bins[b][binsz[b]++] = n;
		}
		return bins;
	}
}
