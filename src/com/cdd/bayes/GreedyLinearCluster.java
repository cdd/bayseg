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

import java.lang.*;
import java.util.*;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.*;
import org.openscience.cdk.fingerprint.model.*;
import org.openscience.cdk.interfaces.IAtomContainer;

/*
 * Greedy linear clustering: a subordinate class for use by CompositeModel, in order to reduce the number of
 * entries to a specific size; the objective is to select a subset that has an even sampling of activity values,
 * as well as being structurally diverse based on the fingerprints.
 */
class GreedyLinearCluster
{
	private List<CompositeModel.Entry> entries;
	private int size;

	// ------------ public methods ------------

	public GreedyLinearCluster(List<CompositeModel.Entry> entries, int size)
	{
		this.entries = entries;
		this.size = size;
	}

	// generates a list of 0-based indices that make up the members of the cluster
	public int[] calculate()
	{
		// sort the incoming indices by value
		final int num = entries.size();
		Integer[] valueOrder = new Integer[num];
		for (int n = 0; n < num; n++) valueOrder[n] = n;
		Arrays.sort(valueOrder, new Comparator<Integer>()
		{
			public int compare(Integer i1, Integer i2)
			{
				final double v1 = entries.get(i1).val, v2 = entries.get(i2).val;
				if (v1 < v2) return -1;
				else if (v1 > v2) return 1;
				else return 0;
			}
		});

		// initiate the opt-in with lowest and highest
		boolean[] mask = new boolean[num];
		mask[valueOrder[0]] = true;
		mask[valueOrder[num - 1]] = true;
		int[] lastIdx = new int[10];
		lastIdx[0] = 0;
		lastIdx[1] = num - 1;
		int lastSz = 2;

		// keep adding more until there's enough
		int npass = size / 5;
		final float invPass = 1.0f / (npass - 1);
		int count = 2;
		while (count < size)
		{
			boolean anything = false;
			for (int n = 0; n < npass && count < size; n++)
			{
				int mid = Math.max(1, (int) Math.round((n - 0.5f) * num * invPass));
				int best = -1;
				float lowest = 0;
				for (int i = mid; (i < mid + 10 || best < 0) && i < num; i++) if (!mask[valueOrder[i]])
				{
					float diff = 0;
					for (int j = 0; j < lastSz; j++)
					{
						diff += tanimoto(entries.get(valueOrder[lastIdx[j]]).fp, entries.get(valueOrder[i]).fp);
					}
					diff /= lastSz;
					if (best < 0 || diff < lowest)
					{
						lowest = diff;
						best = i;
					}
				}
				if (best < 0) continue;

				mask[valueOrder[best]] = true;
				if (lastSz >= lastIdx.length)
				{
					for (int i = 0; i < lastSz - 1; i++) lastIdx[i] = lastIdx[i + 1];
					lastSz--;
				}
				lastIdx[lastSz++] = best;
				count++;
				anything = true;
			}
		}

		// package the results
		int[] retidx = new int[count];
		count = 0;
		for (int n = 0; n < num; n++) if (mask[n]) retidx[count++] = n;
		return retidx;
	}

	// ------------ private methods ------------

	// calculates the Tanimoto coefficient for two lists of hash codes: these are assumed to be sorted and unique, which
	// allows the calculation to be done in O(N) time
	private float tanimoto(int[] fp1, int[] fp2)
	{
		int shared = 0, total = 0;
		final int sz1 = fp1.length, sz2 = fp2.length;
		for (int i1 = 0, i2 = 0; i1 < sz1 || i2 < sz2; total++)
		{
			if (i1 == sz1)
			{
				total += sz2 - i2;
				break;
			}
			if (i2 == sz2)
			{
				total += sz1 - i1;
				break;
			}
			final int v1 = fp1[i1], v2 = fp2[i2];
			if (v1 == v2)
			{
				shared++;
				i1++;
				i2++;
			}
			else if (v1 < v2) i1++;
			else i2++;
		}
		return (float) shared / total;
	}
}
