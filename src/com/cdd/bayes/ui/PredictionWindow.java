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

import javafx.scene.control.*;
import javafx.event.*;
import javafx.geometry.*;
import javafx.stage.*;
import javafx.scene.*;
import javafx.scene.canvas.*;
import javafx.scene.control.*;
import javafx.scene.input.*;
import javafx.scene.layout.*;
import javafx.scene.paint.*;
import javafx.scene.shape.*;
import javafx.application.*;
import javafx.beans.value.*;
import javafx.util.*;

import org.json.*;
import org.controlsfx.control.*;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.exception.CDKException;

/*
 * Computes a bunch of predictions for a set of molecules, and displays them in list form, to be perused.
*/

public class PredictionWindow
{
	// ------------ private data ------------	

    private Stage stage;
	private CompositeModel model;
	private List<CompositeModel.Entry> molecules;

    private VBox content = new VBox();
    private Label labelStatus = new Label();
    
	final int PADDING = 4;
	
    private static final class Prediction
    {
    	IAtomContainer mol;
    	float[] pred;
    	CDKException exception = null;
    	int best;
    	float score;
    	Canvas chart;
    }
    private List<Prediction> predictions = new ArrayList<>(), incoming = new ArrayList<>();

    private Object mutex = new String("!");
   	private boolean shutdown = false;
    
	// ------------ public methods ------------	

	public PredictionWindow(Stage stage, CompositeModel model, List<CompositeModel.Entry> molecules)
	{
		this.stage = stage;
		this.model = model;
		this.molecules = molecules;

        stage.setTitle("Composite Model Predictions");

		ScrollPane scroller = new ScrollPane(content);
		scroller.setHbarPolicy(ScrollPane.ScrollBarPolicy.AS_NEEDED);
		scroller.setVbarPolicy(ScrollPane.ScrollBarPolicy.AS_NEEDED);
		scroller.setFitToWidth(true);
		content.setPadding(new Insets(PADDING));
		content.setSpacing(PADDING);
		content.setPrefWidth(Double.MAX_VALUE);
		
		labelStatus.setText("Beginning prediction...");
		content.getChildren().add(labelStatus);

		BorderPane root = new BorderPane();
		//root.setTop(menuBar);
		root.setCenter(scroller);

		Scene scene = new Scene(root, 500, 800, Color.WHITE);

		stage.setScene(scene);

		stage.setOnCloseRequest(event -> shutdown = true);

		//recreateContent();
		
		new Thread(() -> makePredictions()).start();
 	}

	// ------------ private methods ------------	
	
	// run in main thread: adds any new widgets as necessary
	private void updateContent()
	{
		while (true)
		{
			Prediction p = null;
			synchronized (mutex)
			{
				if (incoming.size() == 0) break;
				p = incoming.remove(0);
			}
			
			HBox hbox = new HBox();
			hbox.setSpacing(PADDING);
			hbox.setAlignment(Pos.CENTER_LEFT);
			
			//String note = p.pred != null ? Arrays.toString(p.pred) : p.exception != null ? p.exception.getMessage() : "?";
			//hbox.getChildren().add(new Label(note));
			
			if (p.chart != null) hbox.getChildren().add(p.chart);

			if (p.pred != null) 
			{
				double[] segments = model.getSegments();
				double min = p.best == 0 ? model.getMinVal() : segments[p.best - 1];
				double max = p.best == segments.length ? model.getMaxVal() : segments[p.best];
				String txt = Util.formatDouble(min, 4) + " .. " + Util.formatDouble(max, 4) + String.format(" (%.1f%%)", 100 * p.score);
				hbox.getChildren().add(new Label(txt));
			}
			
			content.getChildren().add(hbox);
			
			predictions.add(p);
		}
		
		synchronized (mutex)
		{
			if (predictions.size() < molecules.size()) 
			{
				String text = "Predicting " + predictions.size() + " of " + molecules.size();
				labelStatus.setText(text);
			}
			else content.getChildren().remove(labelStatus);
		}
	}
	
	// runs in a background thread, and makes predictions for each molecule, and renders the results nicely
	private void makePredictions()
	{
		List<IAtomContainer> queue = new ArrayList<>();
	
		synchronized (mutex)
		{
			for (CompositeModel.Entry entry : molecules) queue.add(entry.mol);
		}
	
		while (queue.size() > 0)
		{
			IAtomContainer mol = queue.remove(0);

			Prediction p = new Prediction();
			p.mol = mol;
			try {p.pred = model.predictBins(mol);}
			catch (CDKException ex) {p.exception = ex;}
			
			if (p.pred != null)
			{
				p.best = 0;
				for (int n = 1; n < p.pred.length; n++) if (p.pred[n] > p.pred[p.best]) p.best = n;
				
				// score = best * (best / sum of all): the "best" value is a probabiliy (0..1), and if all other probabilities are zero, it can stand as-is; to the extent
				//									   that other options are viable, it decreases proportionately
				p.score = Math.max(0, Math.min(1, p.pred[p.best]));
				if (p.score > 0)
				{
    				float denom = 0;
    				for (float f : p.pred) denom += Math.max(0, Math.min(1, f));
    				p.score *= p.score / denom;
				}
			
				p.chart = renderPredictions(p.pred);
			}
			// !! do the prediction
			// ... and other stuff: rendering
				
			synchronized (mutex)
			{
				incoming.add(p);
		        Platform.runLater(() -> updateContent());
			}
		}
	}
	
	private Canvas renderPredictions(float[] pred)
	{
		final int nbins = pred.length;
		final double barW = Math.max(20, 150 / nbins), barH = 80, edge = 5;
		final double totalW = nbins * barW + 2 * edge, totalH = barH;
	
		Canvas canvas = new Canvas(totalW, totalH);
		GraphicsContext gc = canvas.getGraphicsContext2D();
		
		gc.setFill(Color.WHITE);
		gc.fillRect(0, 0, totalW, totalH);
		
		for (int n = 0; n < nbins; n++)
		{
			float f = Math.max(0, Math.min(1, pred[n]));
			double x = edge + barW * n, w = barW;
			double h = (totalH - 1) * f, y = totalH - 0.5 - h;

			int bg = Util.blendRGB(f, 0xFF0000, 0xFFFF00, 0x00FF00);
			gc.setFill(Util.rgbColor(bg));
			gc.fillRect(x, y, w, h);
			
    		gc.setStroke(Color.BLACK);
    		gc.setLineWidth(1);
			gc.strokeRect(x, y, w, h);
		}
		
		gc.setStroke(Color.BLACK);
		gc.setLineWidth(1);
		//gc.strokeLine(0, totalH - 0.5, totalW, totalH - 0.5);
		gc.strokeRect(0.5, 0.5, totalW - 1, totalH - 1);
		
		return canvas;
	}
}
