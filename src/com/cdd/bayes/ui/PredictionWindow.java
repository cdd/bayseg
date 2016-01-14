/*
 * Bayesian Composite Models
 * 
 * (c) 2015-2016 Collaborative Drug Discovery, Inc.
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
import javafx.scene.image.*;
import javafx.application.*;
import javafx.beans.value.*;
import javafx.util.*;

import org.json.*;
import org.controlsfx.control.*;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.renderer.generators.*;
import org.openscience.cdk.renderer.*;
import org.openscience.cdk.renderer.font.*;
import org.openscience.cdk.renderer.visitor.*;

/*
 * Computes a bunch of predictions for a set of molecules, and displays them in list form, to be perused.
*/

public class PredictionWindow
{
	// ------------ private data ------------	

    private Stage stage;
	private CompositeModel model;
	private List<CompositeModel.Entry> molecules;
	private boolean invertDir;

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
    	Canvas chart, diagram;
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
		
		double minVal = model.getMinVal(), maxVal = model.getMaxVal();
		invertDir = minVal > 0 && maxVal / minVal > 15;
		
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
			
			if (p.diagram != null) hbox.getChildren().add(p.diagram);

			if (p.pred != null) 
			{
				VBox vbox = new VBox();
				vbox.setSpacing(PADDING);
			
				if (p.chart != null) vbox.getChildren().add(p.chart);

				double[] segments = model.getSegments();
				double minVal = p.best == 0 ? model.getMinVal() : segments[p.best - 1];
				double maxVal = p.best == segments.length ? model.getMaxVal() : segments[p.best];

				String txt = Util.formatDouble(invertDir ? maxVal : minVal, 4) + " .. " + Util.formatDouble(invertDir ? minVal : maxVal, 4) + String.format(" (%.1f%%)", 100 * p.score);
				vbox.getChildren().add(new Label(txt));
				
				hbox.getChildren().add(vbox);
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
				p.diagram = renderMolecule(p.mol);
			}
				
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
		final double width = nbins * barW + 2 * edge, height = barH;
	
		Canvas canvas = new Canvas(width, height);
		GraphicsContext gc = canvas.getGraphicsContext2D();
		
		gc.setFill(Color.WHITE);
		gc.fillRect(0, 0, width, height);
		
		for (int n = 0; n < nbins; n++)
		{
			int i = invertDir ? nbins - 1 - n : n;
			float f = Math.max(0, Math.min(1, pred[i]));
			double x = edge + barW * n, w = barW;
			double h = (height - 1) * f, y = height - 0.5 - h;

			int bg = Util.blendRGB(f, 0xFF0000, 0xFFFF00, 0x00FF00);
			gc.setFill(Util.rgbColor(bg));
			gc.fillRect(x, y, w, h);
			
    		gc.setStroke(Color.BLACK);
    		gc.setLineWidth(1);
			gc.strokeRect(x, y, w, h);
		}
		
		gc.setStroke(Color.BLACK);
		gc.setLineWidth(1);
		gc.strokeRect(0.5, 0.5, width - 1, height - 1);
		
		return canvas;
	}
	
	private Canvas renderMolecule(IAtomContainer mol)
	{
		int width = 300, height = 200;
	
		Canvas canvas = new Canvas(width, height);
		GraphicsContext gc = canvas.getGraphicsContext2D();
		
		gc.setFill(Color.WHITE);
		gc.fillRect(0, 0, width, height);

        List<IGenerator<IAtomContainer>> generators = new ArrayList<IGenerator<IAtomContainer>>();
        generators.add(new BasicSceneGenerator());
        generators.add(new BasicAtomGenerator());
        generators.add(new BasicBondGenerator());
        generators.add(new AtomNumberGenerator()); 
        generators.add(new ExtendedAtomGenerator()); 
		
        AtomContainerRenderer renderer=new AtomContainerRenderer(generators,new AWTFontManager());
		renderer.getRenderer2DModel().set(AtomNumberGenerator.WillDrawAtomNumbers.class,Boolean.FALSE);

		// render onto the AWT canvas (yucky but necessary)
		java.awt.image.BufferedImage awtimg = new java.awt.image.BufferedImage(width, height, java.awt.image.BufferedImage.TYPE_INT_ARGB);
		java.awt.Graphics2D g = (java.awt.Graphics2D)awtimg.getGraphics();
		g.setRenderingHint(java.awt.RenderingHints.KEY_ANTIALIASING, java.awt.RenderingHints.VALUE_ANTIALIAS_ON);
		g.setRenderingHint(java.awt.RenderingHints.KEY_STROKE_CONTROL, java.awt.RenderingHints.VALUE_STROKE_PURE);
		
		java.awt.Rectangle box = new java.awt.Rectangle(3, 3, width - 6, height - 6);
        renderer.setup(mol, box);
        renderer.paint(mol, new AWTDrawVisitor(g), box, true);

		WritableImage wimg = new WritableImage(width, height);
		javafx.embed.swing.SwingFXUtils.toFXImage(awtimg, wimg);
		gc.drawImage(wimg, 0, 0);
            
		gc.setStroke(Color.BLACK);
		gc.setLineWidth(1);
		gc.strokeRect(0.5, 0.5, width - 1, height - 1);
		
		return canvas;
	}
}
