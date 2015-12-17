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

/*
 * Renders one or more datasets (using the same X-axis) ready to be displayed.
 */
 
public class RenderDatasets
{
	private List<List<CompositeModel.Entry>> datasets;
	private double[] segments;
	private Canvas canvas = null;
	private GraphicsContext gc = null;
	private int width = 600, height = 200;
	
	private boolean useLog = false; // true if it's something like IC50(conc) and should be calibrated by x = -log(x)
	private double minVal, maxVal; // data range (pre-calibration)
	private double minCal, maxCal; // data range (post-calibration, i.e. maybe log) and with an edge factored in
	private double span, invSpan; // (maxCal-minCal) and 1/span

	// ------------ public methods ------------
	
	public RenderDatasets(List<List<CompositeModel.Entry>> datasets, double[] segments)
	{
		this.datasets = datasets;
		this.segments = segments;
	}
	
	public void draw()
	{
		setupRange();
	
		final int nsets = datasets.size();
		double totalW = width, totalH = height * nsets;
		canvas = new Canvas(totalW, totalH + 15);
		gc = canvas.getGraphicsContext2D();
		Font font = new Font(10);
		
		// background
		gc.setFill(Color.WHITE);
		gc.fillRect(0, 0, totalW, totalH);
		
		for (int n = 0; n < nsets; n++)
		{
			float[] axis = new float[width];
			final float gauss = (float)(0.1 * invSpan * width);
			
			for (CompositeModel.Entry entry : datasets.get(n))
			{
				double val = (calibrated(entry.val) - minCal) * invSpan;
				plotGaussian(axis, (float)(val * width), gauss);
			}
			
			int colFill = n == 0 ? 0x00579C : 0x0E7B4B;
			int colEdge = n == 0 ? 0x4097CC : 0x4EBB8B;
			gc.save();
			drawHistogram(axis, Util.rgbColor(colFill), Util.rgbColor(colEdge), 0, n * height);
			gc.restore();
		}
		
		if (segments != null)
		{
			double deltaY = Math.min(20, height / (2 * segments.length)), midY = 0.5 * deltaY;
			for (double v : segments)
			{
    			double x = (calibrated(v) - minCal) * invSpan * width;
    			gc.save();
    			gc.setLineDashes(1, 3);
    			gc.strokeLine(x, 0, x, height);

    			gc.setFont(font);
    			gc.setFill(Color.BLACK);
    			gc.setTextAlign(TextAlignment.LEFT);
    			String str = Util.formatDouble(v, 4);
    			gc.fillText(str, x + 2, midY);
    			midY += deltaY;
    			gc.restore();
			}
			/*int[] counts = new int[segments.length + 1];
			for (CompositeModel.Entry entry : datasets.get(0))
			{
				int bin = 0;
				for (; bin < segments.length; bin++) if (entry.val < segments[bin]) break;
				counts[bin]++;
			}
			for (int n = 0; n < counts.length; n++)
			{
				double v1 = n == 0 ? minCal : segments[n - 1];
				double v2 = n == segments.length ? maxCal : segments[n];
				double x = (calibrated(0.5 * (v1 + v2)) - minCal) * invSpan * width;
				gc.save();
    			gc.setFont(font);
    			gc.setFill(Color.BLACK);
				gc.setTextAlign(TextAlignment.CENTER);
				String str = String.valueOf(counts[n]);
				gc.fillText(str, x, midY - 0.5 * font.getSize());
				midY += deltaY;
				gc.restore();
			}*/
		}
		
		// frame
		gc.setStroke(Color.BLACK);
		gc.setLineWidth(1);
		gc.setLineJoin(StrokeLineJoin.MITER);
		gc.strokeRect(0.5, 0.5, totalW - 1, totalH - 1);
		for (int n = 0; n < nsets - 1; n++)
		{
			double y = (n + 1) * totalH / nsets;
			gc.strokeLine(0, y, totalW, y);
		}
		
		// axis notes
		for (int n = 0; n < 2; n++)
		{
			double v = n == 0 ? minVal : maxVal;
			double x = (calibrated(v) - minCal) * invSpan * width;
			double y = totalH;
			String str = Util.formatDouble(n == 0 ? minVal : maxVal, 4);
			gc.strokeLine(x, y, x, y + 3);
			gc.setFont(font);
			gc.setFill(Color.BLACK);
			gc.setTextAlign(TextAlignment.CENTER);
			gc.fillText(str, x, y + 12);
		}
	}
	
	public boolean usesLog() {return useLog;}
	public Canvas getCanvas() {return canvas;}
	
	// ------------ private methods ------------

	private void setupRange()
	{
		minVal = Double.POSITIVE_INFINITY;
		maxVal = Double.NEGATIVE_INFINITY;
		for (List<CompositeModel.Entry> dset : datasets) for (CompositeModel.Entry entry : dset)
		{
			minVal = Math.min(minVal, entry.val);
			maxVal = Math.max(maxVal, entry.val);
		}
		
		useLog = minVal > 0 && maxVal / minVal > 15;
		minCal = Math.min(calibrated(minVal), calibrated(maxVal));
		maxCal = Math.max(calibrated(minVal), calibrated(maxVal));
		
		double edge = (maxCal - minCal) * 0.03;
		minCal -= edge;
		maxCal += edge;
		span = maxCal - minCal;
		invSpan = 1 / span;
	}
	
	private double calibrated(double v)
	{
		return useLog ? -Math.log(v) : v;
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

	private void drawHistogram(float[] axis, Color fill, Color edge, double x0, double y0)
	{
		float maxH = axis[0];
		for (int n = 1; n < axis.length; n++) maxH = Math.max(maxH, axis[n]);
		double yscale = (height - 2) / maxH;
	
		for (int pass = 0; pass < 2; pass++)
		{
			gc.beginPath();
    		for (int n = 0; n < axis.length; n++)
    		{
    			double x = x0 + n;
    			double y = y0 + height - axis[n] * yscale - 1;
    			if (n == 0) gc.moveTo(x, y); else gc.lineTo(x, y);
    		}
    		if (pass == 0)
    		{
    			gc.lineTo(x0 + width, y0 + height);
    			gc.lineTo(x0, y0 + height);
    			gc.setFill(fill);
    			gc.fill();
    		}
    		else
    		{
        		gc.setStroke(edge);
        		gc.setLineWidth(1);
        		gc.setLineJoin(StrokeLineJoin.ROUND);
        		gc.stroke();
    		}
		}
	}
}

