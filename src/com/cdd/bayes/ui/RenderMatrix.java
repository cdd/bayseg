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

package com.cdd.bayes.ui;

import com.cdd.bayes.*;
import com.cdd.bayes.util.*;

import java.io.*;
import java.util.*;

import javafx.event.*;
import javafx.geometry.*;
import javafx.stage.*;
import javafx.scene.*;
import javafx.scene.canvas.*;
import javafx.scene.paint.*;
import javafx.scene.shape.*;
import javafx.scene.text.*;

import org.openscience.cdk.exception.CDKException;

/*
 * Renders one or more datasets (using the same X-axis) ready to be displayed.
 */
 
public class RenderMatrix
{
	private CompositeModel model;
	private List<CompositeModel.Entry> dataset;
	private boolean isTesting;
	private boolean invertDir;
	
	private Canvas canvas = null;
	private GraphicsContext gc = null;
	private int sizeWH = 300; // (square)
	
	private double[] segments;
	private int[] binsz; // occupancy for each bin
	private int[][] matrix; // hit/miss matrix
	private int[] offcounts; // hit/miss counts
	private float[] offportion, offrandom, enrichment;

	// ------------ public methods ------------
	
	public RenderMatrix(CompositeModel model, List<CompositeModel.Entry> dataset, boolean isTesting, boolean invertDir)
	{
		this.model = model;
		this.dataset = dataset;
		this.isTesting = isTesting;
		this.invertDir = invertDir;
	}
	
	public void draw()
	{
		final int titleH = 15, footerH = 15;
		canvas = new Canvas(sizeWH, sizeWH + titleH + footerH);
		gc = canvas.getGraphicsContext2D();

		try {setupMatrix();}
		catch (CDKException ex)
		{
			ex.printStackTrace();
			return;
		}
	
		Font font = new Font(10);

		String title = !isTesting ? "Training Matrix" : "Testing Matrix";
		gc.save();
		gc.setFont(font);
		gc.setFill(Color.BLACK);
		gc.setTextAlign(TextAlignment.LEFT);
		gc.fillText(title, 2, font.getSize());
		gc.restore();
		
		drawMatrix(0, titleH, sizeWH);
	}
	
	public Canvas getCanvas() {return canvas;}
	
	public int[] getOffCounts() {return offcounts;}
	public float[] getOffPortion() {return offportion;}
	public float[] getOffRandom() {return offrandom;}
	public float[] getEnrichment() {return enrichment;}
	
	// ------------ private methods ------------

	private void setupMatrix() throws CDKException
	{
		segments = model.getSegments();
		final int nbins = segments.length + 1;
		
		if (!isTesting)
		{
			int[][] bins = model.getAssignedBins();
			binsz = new int[nbins];
			for (int n = 0; n < nbins; n++) binsz[n] = bins[n].length;
			matrix = model.getValidationMatrix();
		}
		else
		{
			binsz = new int[nbins];
			matrix = new int[nbins][];
			for (int n = 0; n < nbins; n++) matrix[n] = new int[nbins];
			
			for (CompositeModel.Entry entry : dataset)
			{
				float[] pred = model.predictBins(entry.mol);
				
    			int want = 0, got = 0;
    			for (int i = 0; i < segments.length; i++) if (entry.val >= segments[i]) want = i + 1;
    			for (int i = 1; i < nbins; i++) if (pred[i] > pred[got]) got = i;
    			binsz[want]++;
    			matrix[want][got]++;
			}
		}
			
		offcounts = new int[nbins];
		for (int i = 0; i < nbins; i++) for (int j = 0; j < nbins; j++)
		{
			offcounts[Math.abs(i - j)] += matrix[i][j];
		}
	}

	private void drawMatrix(double x0, double y0, double wh)
	{
		Font font = new Font(8);

		final int nbins = binsz.length;
		double[] binw = new double[nbins], binp = new double[nbins];
		double binSum = 0;
		for (int n = 0; n < nbins; n++) {binw[n] = Math.sqrt(binsz[n]); binSum += binw[n];}
		double binScale = wh / binSum;
		for (int n = 0; n < nbins; n++) binw[n] *= binScale;
		if (!invertDir)
		{
			binp[0] = 0;
			for (int n = 1; n < nbins; n++) binp[n] = binp[n - 1] + binw[n - 1];
		}
		else
		{
			binp[0] = wh - binw[0];
			for (int n = 1; n < nbins; n++) binp[n] = binp[n - 1] - binw[n];
			
			//for (int n = nbins - 1; n >= 0; n--) binp[n] = binp[n + 1] + binw
		}

		// white background
		gc.setFill(Color.WHITE);		
		gc.fillRect(x0, y0, wh, wh);

		int rgb1 = 0xFFFFFF;
		int rgb2 = !isTesting ? 0x00579C : 0x0E7B4B;
		int rgb3 = !isTesting ? 0x4097CC : 0x4EBB8B;		

		for (int i = 0; i < nbins; i++) for (int j = 0; j < nbins; j++)
		{
			int pop = matrix[i][j];
			double x = x0 + binp[i], y = y0 + binp[j], w = binw[i], h = binw[j];

			if (pop > 0)
			{
    			float ideal = i == j ? binsz[i] : 0.5f * (binsz[i] + binsz[j]);
    			if (ideal == 0) continue;
    			float fract = Math.max(0, Math.min(1, pop / ideal));
    			int rgb = Util.blendRGB(fract, rgb1, rgb2, rgb3);
    			gc.setFill(Util.rgbColor(rgb));
    			gc.fillRect(x, y, w, h);
			}
			
			String label = String.valueOf(pop);
    		gc.save();
    		gc.setFont(font);
    		gc.setFill(Color.BLACK);
    		gc.setTextAlign(TextAlignment.CENTER);
    		gc.fillText(label, x + 0.5 * w, y + 0.5 * (h + font.getSize()), w);
    		gc.restore();
		}
		
		// draw the segmentation values
		gc.save();
		gc.setFont(font);
		gc.setFill(Color.BLACK);
		gc.setTextAlign(TextAlignment.CENTER);
		double segY = y0 + wh + font.getSize() + 3;
		for (int n = 0; n < segments.length; n++)
		{
			double x = x0 + (invertDir ? binp[n] : binp[n + 1]), y = y0 + wh;
			String str = Util.formatDouble(segments[n], 4);
			gc.fillText(str, x, segY);
		}
		gc.setTextAlign(TextAlignment.LEFT);
		gc.fillText(Util.formatDouble(invertDir ? model.getMaxVal() : model.getMinVal(), 4), x0, segY);
		gc.setTextAlign(TextAlignment.RIGHT);
		gc.fillText(Util.formatDouble(invertDir ? model.getMinVal() : model.getMaxVal(), 4), x0 + wh, segY);
		gc.restore();

		// now create the results for tabulation purposes
		offportion = new float[nbins];
		offrandom = new float[nbins];
		enrichment = new float[nbins];

		int totalCount = 0, totalPop = dataset.size();
		float totalPortion = 0, totalRandom = 0;
		for (int n = 0; n < nbins; n++)
		{
			int count = offcounts[n];
			float portion = (float)count / totalPop;
			float random = n == 0 ? 1.0f / nbins : 2.0f * (nbins - n) / (nbins * nbins);
			totalCount += count;
			totalPortion += portion;
			totalRandom += random;

			offportion[n] = totalPortion;
			offrandom[n] = totalRandom;
			enrichment[n] = totalPortion / totalRandom;
		}
		
		// boundary and dividers
		gc.save();
		gc.setStroke(Color.BLACK);
		gc.setLineWidth(1);
		gc.setLineJoin(StrokeLineJoin.MITER);
		gc.strokeRect(x0 + 0.5, y0 + 0.5, sizeWH - 1, sizeWH - 1);
		for (int n = 0; n < nbins - 1; n++)
		{
			double p = binp[invertDir ? n : n + 1];
			gc.strokeLine(x0 + 0.5, y0 + p + 0.5, x0 + wh - 0.5, y0 + p - 0.5);
			gc.strokeLine(x0 + p + 0.5, y0 + 0.5, x0 + p - 0.5, y0 + wh + 3);
		}
		gc.restore();
	}
	
}

